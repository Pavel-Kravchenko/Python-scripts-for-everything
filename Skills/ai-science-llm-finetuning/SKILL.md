---
name: ai-science-llm-finetuning
description: "*Prerequisites: Tier 3 Module 10 (Deep Learning for Biology), PyTorch basics.*"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/01_LLM_Finetuning.ipynb"
---

# Module T5-01: LLM Fine-tuning

*Source: Course notebook `Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/01_LLM_Finetuning.ipynb`*

# Module T5-01: LLM Fine-tuning

**Tier 5 — Modern AI for Science | Module 01**

*Prerequisites: Tier 3 Module 10 (Deep Learning for Biology), PyTorch basics.*

**GPU required** for the training cells. Theory and pattern cells run on CPU. Run on free-tier Colab (T4 GPU).

---

This module covers fine-tuning large language models for domain-specific instruction following using LoRA, quantization, and the `trl` SFTTrainer.

**By the end you will be able to:**
- Explain LoRA: low-rank adapter math and parameter counts
- Apply 4-bit NF4 quantization with bitsandbytes
- Format instruction datasets using chat templates
- Run SFTTrainer to fine-tune a base model
- Generate synthetic instruction data for domain adaptation

**Attribution:** *Patterns inspired by Unsloth AI and Manuel Faysse fine-tuning tutorials. Uses public instruction datasets from HuggingFace Hub.*

## Why this notebook matters

Fine-tuning large language models is the standard way to adapt general-purpose foundation models to domain-specific tasks. In computational biology this means adapting models for biomedical question-answering, literature mining, clinical note summarization, and scientific hypothesis generation. The techniques here — LoRA, 4-bit quantization, and SFTTrainer — are the practical stack used in virtually all LLM fine-tuning projects as of 2024–2025.

## How to work through this notebook

1. Read each section's explanation before running code — the math matters here.
2. The LoRA parameter count cell and quantization config cell run on CPU. Only Section 5 (SFTTrainer) requires a GPU.
3. For GPU cells: run on free-tier Google Colab (T4 is sufficient for 7B with 4-bit quantization).
4. Keep the synthetic `bio_qa` dataset small while learning the pipeline, then scale up.
5. The SFTTrainer block is shown as a printed pattern — copy it into a Colab cell with GPU to execute.

## Common sticking points

- **LoRA rank vs quality tradeoff**: higher rank (r=64) trains more parameters but is not always better — rank r=16 is a strong default. The benefit of higher rank depends on dataset size and task diversity.
- **Quantization and bfloat16**: NF4 quantizes the frozen base model weights; LoRA adapter weights stay in bfloat16 full precision. These are two separate things — do not confuse them.
- **Chat template mismatch**: loading a model with the wrong chat template causes silent garbage outputs, not errors. Always verify the template against the model card before fine-tuning.
- **Overfitting on small datasets**: with fewer than 500 examples, 1–2 epochs is usually enough. Watch eval loss; if it rises while train loss falls, you are overfitting.
- **Padding token**: most base models (Mistral, Llama) do not have a pad token by default. Setting `tokenizer.pad_token = tokenizer.eos_token` is the standard fix.

## 1. LoRA: Low-Rank Adaptation

Full fine-tuning updates all W ≈ 7B parameters — prohibitive without many A100s. LoRA instead trains two small matrices A and B (rank r) that decompose the weight update:

$$\Delta W = B \cdot A, \quad B \in \mathbb{R}^{d \times r}, \; A \in \mathbb{R}^{r \times k}$$

During inference: W_eff = W_pretrained + (α/r) · B·A

**Parameter savings:**
| Model | Full FT params | LoRA (r=16) params | Reduction |
|---|---|---|---|
| 7B (Mistral, Llama) | 7,000,000,000 | ~6,800,000 | 1000× |
| 13B | 13,000,000,000 | ~12,600,000 | 1000× |

**Which modules to target?**
- Minimum: `q_proj`, `v_proj` (query and value in attention)
- Recommended: + `k_proj`, `o_proj`, `gate_proj`, `up_proj`, `down_proj`

```python
def lora_param_count(model_dim, n_heads, n_layers, rank, target_modules):
    """
    Estimate LoRA parameter count.
    model_dim: hidden size (e.g. 4096 for 7B models)
    n_heads: number of attention heads
    n_layers: number of transformer layers
    """
    params_per_layer = 0
    if "q_proj" in target_modules:
        params_per_layer += 2 * rank * model_dim  # A and B matrices
    if "v_proj" in target_modules:
        params_per_layer += 2 * rank * model_dim
    if "k_proj" in target_modules:
        params_per_layer += 2 * rank * model_dim
    if "o_proj" in target_modules:
        params_per_layer += 2 * rank * model_dim

    total = params_per_layer * n_layers
    return total

# Mistral-7B architecture
total_params = 7_000_000_000
for rank in [8, 16, 32, 64]:
    lora_p = lora_param_count(4096, 32, 32, rank,
                               ["q_proj","v_proj","k_proj","o_proj"])
    pct = lora_p / total_params * 100
    print(f"r={rank:3d}: {lora_p:>12,} trainable params ({pct:.3f}% of 7B)")
```

## 2. Quantization: 4-bit NF4

Quantization reduces the numerical precision of model weights from float32 (4 bytes/param) to 4-bit (0.5 bytes/param) — an 8× memory reduction.

**NF4 (NormalFloat 4-bit):** Designed for normally distributed weights (which transformer weights typically are). Values are mapped to 16 fixed points optimally placed for a normal distribution.

| Format | Bits | Memory for 7B | Quality |
|---|---|---|---|
| float32 | 32 | 28 GB | Reference |
| bfloat16 | 16 | 14 GB | ~same as fp32 |
| INT8 | 8 | 7 GB | Small quality loss |
| NF4 | 4 | 3.5 GB | Small quality loss |

**QLoRA:** Combine 4-bit base model (frozen) + LoRA adapters (bf16, trainable). Train only LoRA; inference dequantizes NF4 on the fly.

```python
# Quantization configuration (runs on CPU; actual loading requires GPU)
try:
    from transformers import BitsAndBytesConfig
    import torch

    bnb_config = BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_quant_type="nf4",          # NormalFloat 4-bit
        bnb_4bit_compute_dtype=torch.bfloat16,  # compute in bf16
        bnb_4bit_use_double_quant=True,     # double quantization for extra savings
    )
    print("BitsAndBytesConfig created successfully")
    print(f"  quant_type: {bnb_config.bnb_4bit_quant_type}")
    print(f"  compute_dtype: {bnb_config.bnb_4bit_compute_dtype}")
    print(f"  double_quant: {bnb_config.bnb_4bit_use_double_quant}")
except ImportError:
    print("bitsandbytes not installed — pattern shown below")
    print("""
from transformers import BitsAndBytesConfig
import torch

bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_compute_dtype=torch.bfloat16,
    bnb_4bit_use_double_quant=True,
)

model = AutoModelForCausalLM.from_pretrained(
    "mistralai/Mistral-7B-v0.1",
    quantization_config=bnb_config,
    device_map="auto",
)
""")
```

## 3. Chat Templates

Instruction-tuned models expect inputs in a specific format. The format varies by model family:

**Mistral / Llama-2 Instruct:**
```
[INST] <<SYS>>
{system}
<</SYS>>

{user} [/INST] {assistant}</s>
```

**Alpaca:**
```
### Instruction:
{instruction}

### Input:
{context}

### Response:
{output}
```

**ChatML (OpenAI-compatible):**
```
<|im_start|>system
{system}<|im_end|>
<|im_start|>user
{user}<|im_end|>
<|im_start|>assistant
{assistant}<|im_end|>
```

```python
def format_mistral(system, user, assistant=None):
    prompt = f"[INST] <<SYS>>\n{system}\n<</SYS>>\n\n{user} [/INST]"
    if assistant:
        prompt += f" {assistant}</s>"
    return prompt

def format_alpaca(instruction, input_text="", output=""):
    prompt = f"### Instruction:\n{instruction}\n"
    if input_text:
        prompt += f"\n### Input:\n{input_text}\n"
    prompt += f"\n### Response:\n{output}"
    return prompt

def format_chatml(system, user, assistant=None):
    prompt = f"<|im_start|>system\n{system}<|im_end|>\n"
    prompt += f"<|im_start|>user\n{user}<|im_end|>\n"
    prompt += "<|im_start|>assistant\n"
    if assistant:
        prompt += f"{assistant}<|im_end|>"
    return prompt

# Examples
sys = "You are a helpful bioinformatics assistant."
usr = "What does RPKM stand for and how does it differ from TPM?"
asst = "RPKM (Reads Per Kilobase per Million) normalizes by library size and gene length. TPM also normalizes by both, but divides by the sum of all RPK values, making samples directly comparable."

print("=== Mistral format ===")
print(format_mistral(sys, usr, asst))
print("\n=== Alpaca format ===")
print(format_alpaca(usr, output=asst))
```

## 4. Dataset Preparation

A supervised fine-tuning (SFT) dataset consists of (instruction, response) pairs. Sources:

| Source type | Example | Size needed |
|---|---|---|
| Public instruction sets | Alpaca, OpenHermes | General capability |
| Domain-specific | Curated biomedical Q&A | Domain adaptation |
| Synthetic | Generated by a strong LLM | Data augmentation |

**Golden rules:**
1. Quality > quantity — 500 high-quality examples outperform 10k noisy ones
2. Diversity — cover all task types you want the model to handle
3. Format consistency — use the same template throughout
4. No contamination — hold out test examples before generating training data

```python
import pandas as pd

# Synthetic bioinformatics instruction dataset
bio_qa = [
    {
        "instruction": "Explain the difference between FASTQ and FASTA file formats.",
        "response": "FASTA stores sequences with a header line starting with '>'. FASTQ additionally stores per-base quality scores in Phred format, encoded as ASCII characters. FASTQ is the standard output format for sequencing platforms."
    },
    {
        "instruction": "What is the purpose of PCR duplicate marking in NGS pipelines?",
        "response": "PCR duplicates arise when the same DNA fragment is amplified and sequenced multiple times. Marking them prevents double-counting in variant calling and expression quantification, which would artificially inflate confidence in erroneous calls."
    },
    {
        "instruction": "Describe the three ontologies in Gene Ontology (GO).",
        "response": "Gene Ontology has three orthogonal ontologies: (1) Molecular Function (MF) — the biochemical activity of a gene product; (2) Biological Process (BP) — the pathway or larger process the gene product contributes to; (3) Cellular Component (CC) — where in the cell the gene product is active."
    },
    {
        "instruction": "What is the Bonferroni correction and when should you use it?",
        "response": "The Bonferroni correction divides the significance threshold α by the number of tests m, setting the per-test threshold to α/m. Use it when all tests are independent and you want strict family-wise error rate control. For correlated tests (e.g. genomic variants in LD), the Benjamini-Hochberg FDR procedure is usually preferred."
    },
    {
        "instruction": "Explain what a p-value is and what it is NOT.",
        "response": "A p-value is the probability of observing a test statistic at least as extreme as the one computed, assuming the null hypothesis is true. It is NOT: the probability that H0 is true, the probability that results occurred by chance, or a measure of effect size. Small p-values indicate the data are unlikely under H0, not that the effect is biologically meaningful."
    },
]

# Format as Alpaca
records = []
for qa in bio_qa:
    text = format_alpaca(qa["instruction"], output=qa["response"])
    records.append({"instruction": qa["instruction"], "response": qa["response"], "text": text})

df = pd.DataFrame(records)
print(f"Dataset: {len(df)} examples")
print(f"\nSample formatted text:\n{'='*60}")
print(df["text"].iloc[0])
print(f"{'='*60}")
print(f"\nAverage tokens (approx): {df['text'].str.len().mean() / 4:.0f}")
```

## 5. SFTTrainer Workflow

The `trl` library's `SFTTrainer` wraps HuggingFace `Trainer` with instruction-tuning conveniences: automatic packing, dataset text field selection, and LoRA integration.

```python
# Full SFTTrainer pattern (requires GPU; shown as executable template)
SFTTRAINER_PATTERN = '''
# ━━━━━ Install (Colab) ━━━━━
# !pip install -q unsloth trl peft bitsandbytes transformers accelerate

from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig, TrainingArguments
from peft import LoraConfig, get_peft_model, TaskType
from trl import SFTTrainer
from datasets import Dataset
import torch

MODEL_NAME = "mistralai/Mistral-7B-v0.1"

# 1. Load quantized model
bnb_config = BitsAndBytesConfig(load_in_4bit=True, bnb_4bit_quant_type="nf4",
                                 bnb_4bit_compute_dtype=torch.bfloat16)
model = AutoModelForCausalLM.from_pretrained(MODEL_NAME, quantization_config=bnb_config,
                                              device_map="auto")
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
tokenizer.pad_token = tokenizer.eos_token

# 2. Apply LoRA
lora_cfg = LoraConfig(r=16, lora_alpha=32, target_modules=["q_proj","v_proj","k_proj","o_proj"],
                       lora_dropout=0.05, bias="none", task_type=TaskType.CAUSAL_LM)
model = get_peft_model(model, lora_cfg)
model.print_trainable_parameters()

# 3. Prepare dataset
dataset = Dataset.from_list(records)  # records has "text" field

# 4. Train
training_args = TrainingArguments(
    output_dir="./results", num_train_epochs=2,
    per_device_train_batch_size=2, gradient_accumulation_steps=4,
    learning_rate=2e-4, fp16=True, logging_steps=10,
    gradient_checkpointing=True,
)
trainer = SFTTrainer(model=model, args=training_args, train_dataset=dataset,
                      dataset_text_field="text", max_seq_length=2048, tokenizer=tokenizer)
trainer.train()
trainer.model.save_pretrained("./bio_assistant")
'''

print("SFTTrainer pattern:")
print(SFTTRAINER_PATTERN)
```
