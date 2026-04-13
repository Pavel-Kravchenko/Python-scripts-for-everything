---
name: ai-science-vision-rag
description: "*Prerequisites: Module T5-01 (LLM Fine-tuning), Tier 3 Module 10 (Deep Learning).*"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/02_Vision_RAG/02_Vision_RAG.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pytorch 2.2+, transformers 4.38+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module T5-02: Vision RAG

*Source: Course notebook `Tier_5_Modern_AI_for_Science/02_Vision_RAG/02_Vision_RAG.ipynb`*


**Tier 5 — Modern AI for Science | Module 02**

*Prerequisites: Module T5-01 (LLM Fine-tuning), Tier 3 Module 10 (Deep Learning).*

GPU optional for inference patterns. Theory cells run on CPU.

---

**By the end you will be able to:**
- Explain VLM architecture (encoder + LLM decoder)
- Understand ColPali late-interaction document retrieval
- Implement a page-level RAG pipeline for PDF documents
- Run Qwen2-VL inference patterns
- Evaluate retrieval recall and generation faithfulness

**Attribution:** *Patterns inspired by Unsloth AI and Manuel Faysse Vision RAG tutorials. Uses public PDF documents (arXiv papers).*

## Why this notebook matters

Scientific literature is increasingly stored in PDF format — papers, clinical reports, supplemental materials — where critical content lives in figures, tables, and complex layouts that text extraction breaks. Vision RAG replaces the fragile OCR → text chunking pipeline with a direct image-level retrieval and VLM-based generation pipeline. This is now the standard approach for document-level question answering in research settings.

## How to work through this notebook

1. All cells run on CPU — the retrieval and evaluation sections use numpy simulations.
2. The Qwen2-VL inference pattern (Section 5) requires a GPU to execute; it is shown as a printed template.
3. Work through the MaxSim scoring derivation in Section 3 carefully — it is the core of ColPali.
4. The evaluation metrics in Section 6 are the ones used in the ViDoRe benchmark; implement them in the exercises.

## Common sticking points

- **Single-vector vs late-interaction**: a single embedding per page (like CLIP) loses spatial/layout information. ColPali's per-patch embeddings fix this by letting each query token match independently against any page patch.
- **Why 14×14 patches?** A 224×224 image divided into 16×16 pixel patches gives 196 patches (14×14 grid). Qwen2-VL uses dynamic resolution, so the patch count scales with image size.
- **MaxSim is not symmetric**: the query side sums over all query tokens finding their best patch; this is directional (query→page), not bidirectional.

## 1. Vision-Language Model Architecture

A VLM (Vision-Language Model) combines:
1. **Visual encoder** — processes image patches (ViT: Vision Transformer)
2. **Projection layer** — maps visual tokens to LLM embedding space
3. **LLM decoder** — generates text conditioned on visual + text tokens

**Key VLMs (2024–2025):**
| Model | Size | Strengths |
|---|---|---|
| Qwen2-VL | 2B, 7B, 72B | Multi-image, long context, document understanding |
| LLaVA-1.6 | 7B, 13B, 34B | Strong OCR, open weights |
| InternVL2 | 2B–76B | High benchmark scores |
| PaliGemma | 3B | Compact, versatile |

**Document understanding challenge:** PDF pages have complex layouts (multi-column, tables, figures). Standard text OCR loses layout information. VLMs process the page as an image, preserving layout.

## 2. RAG Pipeline: Retrieval + Generation

**Standard text RAG:**
1. Split document into chunks (sentences/paragraphs)
2. Embed chunks → vector store
3. At query time: embed query, retrieve top-k chunks by similarity
4. Concatenate chunks as context → generate answer

**Vision RAG (ColPali-style):**
1. Render each PDF page as an image
2. Embed pages with a multimodal model → page embedding matrix
3. At query time: embed query text → compute MaxSim against page matrices
4. Retrieve top-k pages → feed as images to VLM
5. VLM generates answer from visual context

**Advantage:** No OCR step — layout, tables, and figures are preserved natively.

```python
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"figure.dpi": 110})

def create_synthetic_page(text_content, page_num, width=595, height=842):
    """Create a synthetic PDF page as a PIL Image."""
    img = Image.new("RGB", (width, height), color="white")
    draw = ImageDraw.Draw(img)

    # Title
    draw.rectangle([40, 40, width-40, 90], fill="#f0f0f0")
    draw.text((50, 55), f"Document Page {page_num}", fill="black")

    # Content
    y = 110
    for line in text_content.split("\n")[:20]:
        draw.text((50, y), line[:80], fill="black")
        y += 25
        if y > height - 60:
            break

    # Simulated figure box
    draw.rectangle([50, height//2, width-50, height//2+200], outline="#888", width=2)
    draw.text((width//2-60, height//2+90), f"[Figure {page_num}]", fill="#888")

    return img

# Create synthetic "document" pages
pages_content = [
    "Introduction to RNA-seq\n\nRNA sequencing is a transcriptomic technique\nthat quantifies gene expression genome-wide.\n\nThe standard pipeline involves:\n1. Library preparation\n2. Sequencing (Illumina)\n3. Quality control (FastQC)\n4. Alignment (STAR, HISAT2)\n5. Quantification (featureCounts)\n6. Differential expression (DESeq2)",
    "Methods: Differential Expression\n\nWe used DESeq2 for differential expression.\nNormalization: median of ratios method.\nStatistical model: negative binomial GLM.\n\nFDR threshold: 0.05 (Benjamini-Hochberg)\nLog2 fold-change threshold: 1.0\n\nTotal genes analyzed: 25,000\nDE genes identified: 1,247",
    "Results: Top DE Genes\n\nUpregulated genes (n=743):\n- MYC: log2FC = 3.2, padj = 1e-45\n- VEGFA: log2FC = 2.8, padj = 2e-38\n- HIF1A: log2FC = 2.1, padj = 5e-29\n\nDownregulated genes (n=504):\n- TP53: log2FC = -2.4, padj = 3e-41\n- RB1: log2FC = -1.9, padj = 7e-35",
    "Discussion\n\nOur findings suggest activation of hypoxic\nresponse pathways as evidenced by HIF1A\nand VEGFA upregulation.\n\nThe downregulation of TP53 and RB1\nis consistent with tumor suppressor\nloss in this cancer subtype.\n\nFuture work: single-cell validation.",
    "Methods: GWAS Analysis\n\nGenotype data from 10,000 samples.\nQuality control: MAF > 0.01, HWE p > 1e-6.\nPopulation stratification corrected by\nincluding top 10 principal components.\n\nAssociation test: logistic regression.\nGenome-wide significance: p < 5e-8.\nClumping: r2 < 0.1, window 250kb.",
]

pages = [create_synthetic_page(content, i+1) for i, content in enumerate(pages_content)]
print(f"Created {len(pages)} synthetic document pages")

# Display page thumbnails
fig, axes = plt.subplots(1, len(pages), figsize=(14, 4))
for ax, page, i in zip(axes, pages, range(len(pages))):
    ax.imshow(np.array(page))
    ax.set_title(f"Page {i+1}", fontsize=9)
    ax.axis("off")
plt.suptitle("Synthetic document pages (simulating PDF rendering)", y=1.02)
plt.tight_layout(); plt.show()
```python

## 3. ColPali: Late-Interaction Retrieval

**Single-vector retrieval (CLIP-style):**
```python
query → q_emb (512d)
page  → p_emb (512d)
score = cosine(q_emb, p_emb)
```python
Problem: compresses complex page into one vector → loses spatial detail.

**Late interaction (ColPali/ColBERT-style):**
```python
query → Q matrix (n_tokens × dim)
page  → P matrix (n_patches × dim)
score = Σ_i max_j(Q[i] · P[j])    ← MaxSim
```python
Each query token finds its best matching page patch independently → preserves spatial detail.

**ColPali training:** A PaliGemma-2 backbone fine-tuned contrastively on (query, relevant page) pairs from document retrieval benchmarks (ViDoRe).

```python
# Simulate ColPali-style embeddings and retrieval
DIM = 128
N_Q_TOKENS = 32
N_PATCHES = 196  # 14×14 patches for 224×224 image

def simulate_page_embedding(page_idx, n_patches=N_PATCHES, dim=DIM, seed=None):
    """Simulate page patch embeddings."""
    r = np.random.default_rng(seed or page_idx)
    emb = r.normal(0, 1, (n_patches, dim))
    emb /= np.linalg.norm(emb, axis=1, keepdims=True)
    return emb

def simulate_query_embedding(query_text, n_tokens=N_Q_TOKENS, dim=DIM):
    """Simulate query token embeddings (text → tokens → embeddings)."""
    seed = sum(ord(c) for c in query_text) % 10000
    r = np.random.default_rng(seed)
    emb = r.normal(0, 1, (n_tokens, dim))
    emb /= np.linalg.norm(emb, axis=1, keepdims=True)
    return emb

def maxsim_score(query_emb, page_emb):
    """ColBERT MaxSim: for each query token, find max similarity to any patch."""
    # query_emb: (n_q, dim), page_emb: (n_p, dim)
    scores = query_emb @ page_emb.T  # (n_q, n_p)
    return scores.max(axis=1).sum()  # sum of per-token maxima

# Simulate relevance: page 2 is relevant to "differential expression" query
page_embeddings = [simulate_page_embedding(i) for i in range(5)]

# Boost page 2 (DESeq2) to be relevant to DE query
relevant_features = simulate_query_embedding("differential expression DESeq2").mean(0)
page_embeddings[1][:10] += relevant_features[:DIM] * 3
page_embeddings[1][:10] /= np.linalg.norm(page_embeddings[1][:10], axis=1, keepdims=True)

# Retrieval
query = "What normalization method was used for differential expression?"
query_emb = simulate_query_embedding(query)
scores = [maxsim_score(query_emb, pe) for pe in page_embeddings]
ranked = sorted(range(5), key=lambda i: -scores[i])

print(f"Query: '{query}'")
print(f"\nRetrieval scores:")
for rank, page_idx in enumerate(ranked):
    print(f"  Rank {rank+1}: Page {page_idx+1} (score = {scores[page_idx]:.3f})")
print(f"\nTop-1 page: {ranked[0]+1} (correct answer: Page 2)")
print(f"Recall@1: {int(ranked[0] == 1)}")
```python

## 4. Full RAG Pipeline

```python
def simple_rag_pipeline(query, pages, page_embeddings, top_k=2, verbose=True):
    """
    Simplified RAG pipeline:
    1. Embed query
    2. Retrieve top-k pages by MaxSim
    3. (Simulate) VLM generation from retrieved pages
    """
    # Step 1: embed query
    query_emb = simulate_query_embedding(query)

    # Step 2: retrieve
    scores = [maxsim_score(query_emb, pe) for pe in page_embeddings]
    top_idx = sorted(range(len(scores)), key=lambda i: -scores[i])[:top_k]

    if verbose:
        print(f"Query: '{query}'")
        print(f"Retrieved pages: {[i+1 for i in top_idx]}")
        print(f"Scores: {[f'{scores[i]:.2f}' for i in top_idx]}")

    # Step 3: simulate VLM generation
    # In practice: feed pages[top_idx] as images to Qwen2-VL with the query
    page_contents_local = [pages_content[i] for i in top_idx]
    context_summary = " | ".join([c[:80] for c in page_contents_local])

    return {
        "retrieved_pages": [i+1 for i in top_idx],
        "scores": [scores[i] for i in top_idx],
        "simulated_context": context_summary,
    }

# Test queries
queries = [
    "What normalization method was used?",
    "Which genes were upregulated?",
    "What was the GWAS sample size?",
]

for q in queries:
    result = simple_rag_pipeline(q, pages, page_embeddings, top_k=2)
    print()
```python

## 5. Qwen2-VL Inference Pattern

Qwen2-VL supports multiple images in a single prompt, making it ideal for multi-page document understanding. It uses a dynamic resolution approach — patches scale with image size.

```python
QWEN_PATTERN = """
# ━━━ Install (Colab) ━━━
# !pip install -q transformers accelerate qwen-vl-utils

from transformers import Qwen2VLForConditionalGeneration, AutoProcessor
from qwen_vl_utils import process_vision_info
import torch

model = Qwen2VLForConditionalGeneration.from_pretrained(
    "Qwen/Qwen2-VL-7B-Instruct",
    torch_dtype=torch.bfloat16,
    device_map="auto",
    attn_implementation="flash_attention_2",
)
processor = AutoProcessor.from_pretrained(
    "Qwen/Qwen2-VL-7B-Instruct",
    min_pixels=256*28*28,
    max_pixels=1280*28*28,
)

def answer_from_pages(pages, question):
    \"\"\"Answer a question using retrieved document pages.\"\"\"
    content = [{"type": "image", "image": page} for page in pages]
    content.append({"type": "text", "text": question})

    messages = [{"role": "user", "content": content}]
    text = processor.apply_chat_template(messages, tokenize=False,
                                          add_generation_prompt=True)
    image_inputs, video_inputs = process_vision_info(messages)
    inputs = processor(text=[text], images=image_inputs, padding=True,
                       return_tensors="pt").to(model.device)

    with torch.no_grad():
        ids = model.generate(**inputs, max_new_tokens=512, temperature=0.1, do_sample=False)

    # Decode only the new tokens
    answer = processor.decode(ids[0][inputs.input_ids.shape[1]:], skip_special_tokens=True)
    return answer
"""

print("Qwen2-VL inference pattern:")
print(QWEN_PATTERN)
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
