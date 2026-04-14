---
name: ai-science-llm-finetuning
description: LLM fine-tuning with HuggingFace Transformers — LoRA math, QLoRA 4-bit NF4 quantization, chat templates, SFTTrainer workflow.
tool_type: python
primary_tool: HuggingFace Transformers
---

# LLM Fine-tuning (LoRA / QLoRA)

GPU required for training. LoRA parameter count and quantization config cells run on CPU.

## Pitfalls

- **LoRA rank vs quality**: Higher rank (r=64) trains more parameters but is not always better — r=16 is the strong default. Benefit of higher rank depends on dataset size and task diversity.
- **Quantization ≠ adapter precision**: NF4 quantizes the frozen base model weights; LoRA adapter weights stay in bfloat16. These are separate concerns.
- **Chat template mismatch**: Loading a model with the wrong chat template causes silent garbage outputs, not errors. Always verify template against the model card.
- **Overfitting on small datasets**: <500 examples → 1–2 epochs. Watch eval loss; if it rises while train loss falls, stop early.
- **Missing pad token**: Most base models (Mistral, Llama) have no pad token. Set `tokenizer.pad_token = tokenizer.eos_token`.
- **Quality > quantity**: 500 high-quality examples outperform 10k noisy ones. Diversity and format consistency matter more than volume.

## LoRA: Parameter Math

$$\Delta W = B \cdot A, \quad B \in \mathbb{R}^{d \times r}, \; A \in \mathbb{R}^{r \times k}$$

During inference: `W_eff = W_pretrained + (α/r) · B·A`

| Model | Full FT params | LoRA r=16 | Reduction |
|---|---|---|---|
| 7B (Mistral/Llama) | 7,000,000,000 | ~6,800,000 | ~1000× |
| 13B | 13,000,000,000 | ~12,600,000 | ~1000× |

**Modules to target:**
- Minimum: `q_proj`, `v_proj`
- Recommended: + `k_proj`, `o_proj`, `gate_proj`, `up_proj`, `down_proj`

## Memory by Precision

| Format | Bits | 7B model | Quality |
|---|---|---|---|
| float32 | 32 | 28 GB | Reference |
| bfloat16 | 16 | 14 GB | ~same |
| INT8 | 8 | 7 GB | Small loss |
| NF4 | 4 | 3.5 GB | Small loss |

## Chat Templates

```python
# Mistral / Llama-2 Instruct
"[INST] <<SYS>>\n{system}\n<</SYS>>\n\n{user} [/INST] {assistant}</s>"

# Alpaca
"### Instruction:\n{instruction}\n\n### Input:\n{context}\n\n### Response:\n{output}"

# ChatML (OpenAI-compatible)
"<|im_start|>system\n{system}<|im_end|>\n<|im_start|>user\n{user}<|im_end|>\n<|im_start|>assistant\n{assistant}<|im_end|>"
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
    return prompt + f"\n### Response:\n{output}"
```

## SFTTrainer Workflow

```python
# pip install unsloth trl peft bitsandbytes transformers accelerate
from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig, TrainingArguments
from peft import LoraConfig, get_peft_model, TaskType
from trl import SFTTrainer
from datasets import Dataset
import torch

MODEL_NAME = "mistralai/Mistral-7B-v0.1"

# 1. Load quantized base model
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_compute_dtype=torch.bfloat16,
    bnb_4bit_use_double_quant=True,  # extra memory savings
)
model = AutoModelForCausalLM.from_pretrained(
    MODEL_NAME, quantization_config=bnb_config, device_map="auto"
)
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
tokenizer.pad_token = tokenizer.eos_token  # required for most base models

# 2. Apply LoRA
lora_cfg = LoraConfig(
    r=16, lora_alpha=32,
    target_modules=["q_proj", "v_proj", "k_proj", "o_proj"],
    lora_dropout=0.05, bias="none", task_type=TaskType.CAUSAL_LM
)
model = get_peft_model(model, lora_cfg)
model.print_trainable_parameters()

# 3. Prepare dataset — each record needs a "text" field with formatted prompt
dataset = Dataset.from_list(records)

# 4. Train
training_args = TrainingArguments(
    output_dir="./results", num_train_epochs=2,
    per_device_train_batch_size=2, gradient_accumulation_steps=4,
    learning_rate=2e-4, fp16=True, logging_steps=10,
    gradient_checkpointing=True,
)
trainer = SFTTrainer(
    model=model, args=training_args, train_dataset=dataset,
    dataset_text_field="text", max_seq_length=2048, tokenizer=tokenizer
)
trainer.train()
trainer.model.save_pretrained("./bio_assistant")
```

## Dataset Quality Rules

| Rule | Detail |
|---|---|
| Quality > quantity | 500 curated > 10k noisy |
| Diversity | Cover all task types |
| Format consistency | One template throughout |
| No contamination | Hold out test set before generating training data |
