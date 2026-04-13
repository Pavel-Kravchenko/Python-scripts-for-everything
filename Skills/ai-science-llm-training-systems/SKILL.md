---
name: ai-science-llm-training-systems
description: "**Tier 5 — Modern AI for Science | Module 01 · Notebook 2**"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/02_llm_training_systems.ipynb"
primary_tool: Pandas
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.1+, pytorch 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module T5-01B: LLM Training Systems (Tracking, Epochs, and Ablations)

*Source: Course notebook `Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/02_llm_training_systems.ipynb`*


**Tier 5 — Modern AI for Science | Module 01 · Notebook 2**

*Prerequisites: 01_LLM_Finetuning.ipynb*

---

**By the end of this notebook you will be able to:**
1. Set up an experiment system for LLM fine-tuning
2. Track training/validation metrics per epoch and step
3. Design and analyze ablation studies
4. Build a reproducible run registry with configs and results
5. Create a practical decision loop for choosing the next experiment

## Why this notebook matters

Reproducibility in LLM fine-tuning is harder than in classical ML — runs are expensive, non-deterministic at scale, and results depend on subtle interactions between dataset, template, learning rate schedule, and rank. A proper experiment system (config tracking, ablation matrix, run registry) is what separates a one-off fine-tune from a repeatable research or production pipeline.

## How to work through this notebook

1. All cells are CPU-friendly — this notebook is about infrastructure, not training.
2. Run sections in order; each section builds on the `TrainConfig` dataclass from Section 1.
3. The tracking backend templates (W&B, MLflow, TensorBoard) are shown as code comments — uncomment and adapt for your environment.
4. Focus on the ablation study design pattern in Section 4 before trying to optimize any hyperparameter.

## Common sticking points

- **Why simulated training curves?** Real training requires GPU hours. The simulated curves here let you practice the analysis workflow (early stopping, epoch selection, ablation comparison) without spending compute.
- **One factor at a time**: changing learning rate and rank simultaneously makes it impossible to attribute differences. Design clean ablations.
- **Random seed discipline**: fix seeds for dataset splits, weight initialization, and data shuffling before starting a run. Adding seeds retroactively breaks reproducibility.
- **Primary vs secondary metrics**: pick one primary metric (e.g., eval_loss) before running. Cherry-picking secondary metrics post-hoc is p-hacking for LLMs.

## 1. What an LLM training system includes

A robust setup has five parts:
1. **Config management** (all hyperparameters captured)
2. **Tracking backend** (W&B / MLflow / TensorBoard)
3. **Evaluation protocol** (fixed validation + task metrics)
4. **Ablation framework** (controlled experiments)
5. **Run registry** (what changed, why, and what improved)

Without this, you cannot reliably tell whether improvement came from model changes or noise.

```python
from dataclasses import dataclass, asdict
import pandas as pd
import numpy as np

np.random.seed(42)

@dataclass
class TrainConfig:
    run_name: str
    model_name: str
    lora_rank: int
    learning_rate: float
    batch_size: int
    epochs: int
    warmup_ratio: float
    weight_decay: float
    max_seq_len: int
    seed: int

cfg = TrainConfig(
    run_name='baseline-lora-r16',
    model_name='mistral-7b-instruct',
    lora_rank=16,
    learning_rate=2e-4,
    batch_size=8,
    epochs=3,
    warmup_ratio=0.03,
    weight_decay=0.01,
    max_seq_len=2048,
    seed=42,
)

pd.Series(asdict(cfg))
```python

## 2. Tracking setup templates

Use one tracker consistently across all runs in a project.

### Weights & Biases template
```python
# import wandb
# wandb.init(project='llm-finetuning', name=cfg.run_name, config=asdict(cfg))
# wandb.log({'train/loss': loss, 'eval/loss': eval_loss, 'eval/accuracy': acc}, step=global_step)
# wandb.finish()
```python

### MLflow template
```python
# import mlflow
# mlflow.set_experiment('llm-finetuning')
# with mlflow.start_run(run_name=cfg.run_name):
#     mlflow.log_params(asdict(cfg))
#     mlflow.log_metric('train_loss', float(loss), step=global_step)
#     mlflow.log_metric('eval_loss', float(eval_loss), step=global_step)
```python

### TensorBoard template
```python
# from torch.utils.tensorboard import SummaryWriter
# writer = SummaryWriter(log_dir=f'runs/{cfg.run_name}')
# writer.add_scalar('train/loss', loss, global_step)
# writer.add_scalar('eval/loss', eval_loss, global_step)
# writer.close()
```python

## 3. Epoch-level tracking and early stopping logic

```python
def simulate_training(epochs=5, base_train=2.2, base_eval=2.4, noise=0.03):
    rows = []
    train_loss = base_train
    eval_loss = base_eval
    for ep in range(1, epochs + 1):
        train_loss = train_loss * (0.86 + np.random.uniform(-0.02, 0.02))
        eval_loss = eval_loss * (0.90 + np.random.uniform(-0.03, 0.03)) + np.random.normal(0, noise)
        eval_acc = max(0.0, min(1.0, 1.35 - eval_loss / 2.0 + np.random.normal(0, 0.01)))
        rows.append({'epoch': ep, 'train_loss': train_loss, 'eval_loss': eval_loss, 'eval_acc': eval_acc})
    return pd.DataFrame(rows)

history = simulate_training(epochs=8)
history
```python

```python
def best_epoch(df: pd.DataFrame):
    i = df['eval_loss'].idxmin()
    return int(df.loc[i, 'epoch']), float(df.loc[i, 'eval_loss']), float(df.loc[i, 'eval_acc'])

ep, loss, acc = best_epoch(history)
print(f'Best epoch = {ep}, eval_loss = {loss:.4f}, eval_acc = {acc:.4f}')
```python

## 4. Ablation study design

Good ablations change **one factor at a time**.

Example factors:
- LoRA rank: 8 vs 16 vs 32
- Learning rate: 1e-4 vs 2e-4 vs 5e-4
- Sequence length: 1024 vs 2048
- Warmup ratio: 0.01 vs 0.03 vs 0.1

```python
def simulate_run(lora_rank, lr, seq_len):
    # toy score surface: higher rank helps a bit, too-high lr hurts, longer seq helps slightly
    score = 0.64
    score += 0.015 * np.log2(lora_rank)
    score -= 0.10 * abs(np.log10(lr) - np.log10(2e-4))
    score += 0.01 if seq_len == 2048 else 0.0
    score += np.random.normal(0, 0.004)
    eval_loss = 1.4 - score + np.random.normal(0, 0.01)
    return score, eval_loss

rows = []
for r in [8, 16, 32]:
    for lr in [1e-4, 2e-4, 5e-4]:
        for s in [1024, 2048]:
            score, ev = simulate_run(r, lr, s)
            rows.append({'lora_rank': r, 'lr': lr, 'seq_len': s, 'eval_score': score, 'eval_loss': ev})

ablation_df = pd.DataFrame(rows).sort_values('eval_score', ascending=False)
ablation_df.head(10)
```python

```python
summary = ablation_df.groupby('lora_rank', as_index=False)['eval_score'].mean().sort_values('eval_score', ascending=False)
summary
```python

## 5. Run registry and experiment memory

Save minimal metadata for every run:
- run_id, date, git commit, dataset version
- config values
- best metrics
- note: what changed vs previous run

```python
registry = pd.DataFrame([
    {'run_id': 'r001', 'commit': 'abc123', 'change': 'baseline r16', 'best_eval_loss': 0.694, 'best_eval_acc': 0.712},
    {'run_id': 'r002', 'commit': 'def456', 'change': 'rank 32',     'best_eval_loss': 0.681, 'best_eval_acc': 0.724},
    {'run_id': 'r003', 'commit': 'ghi789', 'change': 'lr 5e-4',     'best_eval_loss': 0.741, 'best_eval_acc': 0.689},
])
registry.sort_values('best_eval_loss')
```python

## 6. Practical checklist for real training

Before launching:
1. Freeze data splits and seed
2. Log full config + package versions
3. Define primary metric (e.g., eval_loss) and secondary metric (task F1/accuracy)
4. Set checkpoint cadence and early-stop policy
5. Predefine ablation matrix and stopping budget

After run:
1. Compare against best prior run
2. Record interpretation (why change helped/hurt)
3. Decide next run from evidence, not intuition alone

## 7. Suggested stack for big-LLM training

- **Trainer/runtime**: Hugging Face `transformers` + `trl` + `accelerate`
- **Scaling**: DeepSpeed or FSDP for multi-GPU/multi-node
- **Tracking**: W&B or MLflow
- **Data/versioning**: DVC or immutable dataset snapshots
- **Orchestration**: Slurm/Kubernetes + reproducible launch scripts

Start with one-node reproducible runs, then scale out.

## Validated Resources

- [Hugging Face Transformers](https://github.com/huggingface/transformers)
- [TRL (SFTTrainer, DPO, etc.)](https://github.com/huggingface/trl)
- [Accelerate](https://github.com/huggingface/accelerate)
- [DeepSpeed](https://github.com/microsoft/DeepSpeed)
- [PyTorch FSDP docs](https://pytorch.org/docs/stable/fsdp.html)
- [Weights & Biases docs](https://docs.wandb.ai/)
- [MLflow docs](https://mlflow.org/docs/latest/index.html)
- [TensorBoard docs](https://www.tensorflow.org/tensorboard)

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
