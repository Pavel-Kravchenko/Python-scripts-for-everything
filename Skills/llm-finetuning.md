---
name: llm-finetuning
description: LLM fine-tuning — LoRA adapter math, 4-bit NF4 quantization with bitsandbytes, chat template formatting, SFTTrainer workflow with trl, synthetic instruction data generation
---

## When to Use

Use this skill when:
- Fine-tuning a pre-trained LLM on a custom instruction/chat dataset
- Applying LoRA (Low-Rank Adaptation) to reduce trainable parameters
- Quantizing a model to 4-bit or 8-bit for memory efficiency
- Designing chat templates (system/user/assistant format)
- Building synthetic instruction datasets for domain adaptation
- Evaluating fine-tuned models on held-out prompts

## Quick Reference

| Concept | Value / Pattern | Notes |
|---|---|---|
| LoRA rank r | 8–64 | Higher r = more parameters, better fit, more memory |
| LoRA alpha | 2× rank typical | Scaling factor: effective_lr = alpha/r × lr |
| Target modules | `q_proj, v_proj` (min) | Add `k_proj, o_proj, gate_proj` for more capacity |
| Quantization | NF4 (4-bit) or INT8 | NF4 best for quality; INT8 faster inference |
| Max seq length | 2048–4096 | Longer = more memory; use `max_seq_length` in SFTTrainer |
| Batch size | 1–4 with grad accum 4–16 | Effective batch = batch × grad_accum |
| Learning rate | 2e-4 to 2e-5 | Cosine schedule with warmup |
| Epochs | 1–3 | Overfitting risk high with small datasets |

## Key Patterns

**Pattern 1: Load quantized base model**
```python
from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig
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
tokenizer = AutoTokenizer.from_pretrained("mistralai/Mistral-7B-v0.1")
tokenizer.pad_token = tokenizer.eos_token
```

**Pattern 2: Apply LoRA with PEFT**
```python
from peft import LoraConfig, get_peft_model, TaskType

lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    target_modules=["q_proj", "v_proj", "k_proj", "o_proj"],
    lora_dropout=0.05,
    bias="none",
    task_type=TaskType.CAUSAL_LM,
)
model = get_peft_model(model, lora_config)
model.print_trainable_parameters()
# "trainable params: 6,815,744 || all params: 3,759,439,872 || trainable%: 0.18"
```

**Pattern 3: Chat template formatting**
```python
def format_instruction(system, user, assistant=None):
    """Mistral/Llama instruction format."""
    prompt = f"[INST] <<SYS>>\n{system}\n<</SYS>>\n\n{user} [/INST]"
    if assistant:
        prompt += f" {assistant}</s>"
    return prompt

# Alpaca format
def format_alpaca(instruction, input_text="", output=""):
    prompt = f"### Instruction:\n{instruction}\n"
    if input_text:
        prompt += f"\n### Input:\n{input_text}\n"
    prompt += f"\n### Response:\n{output}"
    return prompt
```

**Pattern 4: SFTTrainer setup**
```python
from trl import SFTTrainer
from transformers import TrainingArguments

training_args = TrainingArguments(
    output_dir="./results",
    num_train_epochs=2,
    per_device_train_batch_size=2,
    gradient_accumulation_steps=4,
    learning_rate=2e-4,
    lr_scheduler_type="cosine",
    warmup_ratio=0.05,
    fp16=True,
    logging_steps=10,
    save_strategy="epoch",
)
trainer = SFTTrainer(
    model=model,
    args=training_args,
    train_dataset=dataset,
    dataset_text_field="text",
    max_seq_length=2048,
    tokenizer=tokenizer,
)
trainer.train()
```

**Pattern 5: Synthetic data generation**
```python
def generate_qa_pairs(context, model, tokenizer, n=10):
    """Generate synthetic instruction-response pairs from a context."""
    prompt = f"Generate {n} question-answer pairs about:\n{context}\nFormat: Q: ... A: ..."
    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
    outputs = model.generate(**inputs, max_new_tokens=500, temperature=0.7)
    return tokenizer.decode(outputs[0], skip_special_tokens=True)
```

## Code Templates

**Template 1: End-to-end fine-tuning pipeline**
```python
from datasets import Dataset
import pandas as pd

# Prepare dataset
records = [{"text": format_alpaca(instr, inp, out)}
           for instr, inp, out in zip(instructions, inputs, outputs)]
dataset = Dataset.from_list(records)
dataset = dataset.train_test_split(test_size=0.1)

# Fine-tune
trainer = SFTTrainer(
    model=model, args=training_args,
    train_dataset=dataset["train"],
    eval_dataset=dataset["test"],
    dataset_text_field="text",
    max_seq_length=2048,
)
trainer.train()
trainer.model.save_pretrained("finetuned_model")
```

## Common Pitfalls

- **OOM on A100:** reduce batch size + increase gradient accumulation; enable `gradient_checkpointing=True`
- **Tokenizer padding:** always set `pad_token = eos_token` for decoder-only models; left-padding for generation
- **LoRA rank too high:** r > 64 rarely helps and wastes memory; start with r=16
- **Dataset format:** SFTTrainer expects `dataset_text_field` to contain the full formatted text including response; tokenization is automatic
- **Evaluation:** ROUGE/BLEU are poor metrics for instruction following; use LLM-as-judge or task-specific benchmarks

## Related Skills

- `ml-deep-learning-bio` — PyTorch fundamentals and model architecture
- `vision-rag` — extends LLMs with visual document retrieval
- `diffusion-generative` — separate generative paradigm (diffusion vs autoregressive)
