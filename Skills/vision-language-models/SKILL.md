---
name: vision-language-models
description: Vision-language model inference patterns for scientific documents.
primary_tool: HuggingFace Transformers
---

## When to Use

Use this atomic skill for focused work on **vision-language-models** without bundling unrelated topics.

## Quick Reference

This skill was split from `vision-rag.md` to keep topics independent and self-contained.

## Core Patterns

Use the parent material below as the source reference, then keep implementations specific to this topic.

## Source Reference (from merged skill)

name: vision-rag
description: Vision-language models (VLMs), ColPali late-interaction retrieval for documents, RAG pipeline construction, Qwen2-VL inference patterns

## When to Use

Use this skill when:
- Building a retrieval-augmented generation system over PDF/image documents
- Using vision-language models for document understanding
- Implementing ColPali-style late-interaction retrieval
- Running Qwen2-VL or similar VLMs for visual Q&A
- Evaluating retrieval recall and generation faithfulness

## Quick Reference

| Component | Tool / Model | Notes |
|---|---|---|
| Document embedding | ColPali (`vidore/colpali-v1.2`) | Page-level embeddings; late interaction |
| VLM inference | Qwen2-VL, LLaVA, InternVL | Multi-image support in Qwen2-VL |
| PDF rendering | `pdf2image` / `pymupdf` | Render each page as PIL Image |
| Vector store | FAISS, ChromaDB | Store page embeddings |
| Late interaction | MaxSim scoring | ColBERT-style token-level matching |
| Retrieval eval | Recall@k | k=1,5,10 typical |
| Generation eval | Faithfulness (NLI-based) | `TRUE` model or LLM-as-judge |

## Key Patterns

**Pattern 1: Render PDF pages**
```python
from pdf2image import convert_from_path

pages = convert_from_path("paper.pdf", dpi=150)  # returns list of PIL Images
print(f"{len(pages)} pages, size: {pages[0].size}")
```python

**Pattern 2: ColPali embedding**
```python
from transformers import AutoProcessor
from colpali_engine.models import ColPali
import torch

model = ColPali.from_pretrained("vidore/colpali-v1.2", torch_dtype=torch.bfloat16)
processor = AutoProcessor.from_pretrained("vidore/colpali-v1.2")

# Embed pages
def embed_pages(pages):
    inputs = processor(images=pages, return_tensors="pt")
    with torch.no_grad():
        embeddings = model(**inputs)  # (n_pages, n_patches, dim)
    return embeddings

# Embed query
def embed_query(text):
    inputs = processor(text=[text], return_tensors="pt")
    with torch.no_grad():
        return model(**inputs)  # (1, n_tokens, dim)
```python

**Pattern 3: MaxSim retrieval (ColBERT-style)**
```python
def maxsim_score(query_emb, page_emb):
    """Late interaction: max similarity per query token, then sum."""
    scores = torch.einsum("qd,pd->qp", query_emb, page_emb)  # (n_q, n_p)
    return scores.max(dim=1).values.sum().item()

def retrieve(query_emb, page_embeddings, top_k=5):
    scores = [maxsim_score(query_emb[0], pe) for pe in page_embeddings]
    top_idx = sorted(range(len(scores)), key=lambda i: -scores[i])[:top_k]
    return top_idx, [scores[i] for i in top_idx]
```python

**Pattern 4: Qwen2-VL inference**
```python
from transformers import Qwen2VLForConditionalGeneration, AutoProcessor
from qwen_vl_utils import process_vision_info

model = Qwen2VLForConditionalGeneration.from_pretrained(
    "Qwen/Qwen2-VL-7B-Instruct", torch_dtype="auto", device_map="auto"
)
processor = AutoProcessor.from_pretrained("Qwen/Qwen2-VL-7B-Instruct")

messages = [{"role": "user", "content": [
    {"type": "image", "image": page_image},
    {"type": "text",  "text": "What are the main findings reported on this page?"}
]}]
text = processor.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
image_inputs, _ = process_vision_info(messages)
inputs = processor(text=[text], images=image_inputs, return_tensors="pt")
outputs = model.generate(**inputs, max_new_tokens=256)
answer = processor.decode(outputs[0], skip_special_tokens=True)
```python

## Code Templates

**Template 1: Full RAG pipeline**
```python
# 1. Index
pages = convert_from_path("document.pdf", dpi=150)
page_embeddings = [embed_pages([p]) for p in pages]

# 2. Retrieve
query = "What is the main contribution?"
query_emb = embed_query(query)
top_pages, scores = retrieve(query_emb, page_embeddings, top_k=3)

# 3. Generate
context_pages = [pages[i] for i in top_pages]
answer = qwen_answer(context_pages, query)
```python

## Pitfalls

- **Memory:** Qwen2-VL-7B requires ~18 GB GPU RAM in 4-bit; use smaller variant (2B) for CPU
- **DPI tradeoff:** higher DPI = better OCR quality but slower embedding; 150 DPI is a good default
- **Late interaction vs single vector:** ColPali outperforms CLIP-style single-vector retrieval on document pages; use ColPali for multi-column layouts
- **Hallucination:** VLMs hallucinate on low-quality scans; add confidence thresholding

## Related Skills

- `llm-finetuning` — fine-tune the generation component for domain adaptation
- `ml-deep-learning-bio` — transformer architecture foundations


## Related Skills

- `vision-language-models` (this file)
- `vision-rag` (legacy merged skill)
