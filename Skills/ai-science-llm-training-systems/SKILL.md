---
name: ai-science-llm-training-systems
description: "Module T5-01B: LLM Training Systems (Tracking, Epochs, and Ablations) with Pandas"
tool_type: python
primary_tool: Pandas
---

# LLM Training Systems: Tracking, Epochs, and Ablations

## Robust Training Setup

Five required components:
1. **Config management** — all hyperparameters captured
2. **Tracking backend** — W&B / MLflow / TensorBoard
3. **Evaluation protocol** — fixed validation + task metrics
4. **Ablation framework** — controlled experiments
5. **Run registry** — what changed, why, and what improved

**Key discipline rules:**
- Change one factor at a time — simultaneous changes make attribution impossible
- Fix seeds for dataset splits, weight init, and data shuffling before starting any run
- Pick one primary metric (e.g., `eval_loss`) before running — cherry-picking post-hoc is p-hacking

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
```

## Epoch Tracking and Early Stopping

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
```

```python
def best_epoch(df: pd.DataFrame):
    i = df['eval_loss'].idxmin()
    return int(df.loc[i, 'epoch']), float(df.loc[i, 'eval_loss']), float(df.loc[i, 'eval_acc'])

ep, loss, acc = best_epoch(history)
print(f'Best epoch = {ep}, eval_loss = {loss:.4f}, eval_acc = {acc:.4f}')
```

## Ablation Study Design

Ablation factors (one at a time):
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
```

```python
summary = ablation_df.groupby('lora_rank', as_index=False)['eval_score'].mean().sort_values('eval_score', ascending=False)
summary
```

## Run Registry

Minimum fields per run: run_id, date, git commit, dataset version, config values, best metrics, note on what changed.

```python
registry = pd.DataFrame([
    {'run_id': 'r001', 'commit': 'abc123', 'change': 'baseline r16', 'best_eval_loss': 0.694, 'best_eval_acc': 0.712},
    {'run_id': 'r002', 'commit': 'def456', 'change': 'rank 32',     'best_eval_loss': 0.681, 'best_eval_acc': 0.724},
    {'run_id': 'r003', 'commit': 'ghi789', 'change': 'lr 5e-4',     'best_eval_loss': 0.741, 'best_eval_acc': 0.689},
])
registry.sort_values('best_eval_loss')
```

## Pre/Post-Run Checklist

**Before launching:**
1. Freeze data splits and seed
2. Log full config + package versions
3. Define primary metric (e.g., `eval_loss`) and secondary metric (task F1/accuracy)
4. Set checkpoint cadence and early-stop policy
5. Predefine ablation matrix and stopping budget

**After run:**
1. Compare against best prior run
2. Record interpretation (why change helped/hurt)
3. Decide next run from evidence, not intuition alone

## Recommended Stack

- **Trainer/runtime**: HuggingFace `transformers` + `trl` + `accelerate`
- **Scaling**: DeepSpeed or FSDP for multi-GPU/multi-node
- **Tracking**: W&B or MLflow
- **Data/versioning**: DVC or immutable dataset snapshots
- **Orchestration**: Slurm/Kubernetes + reproducible launch scripts

## References

- [Hugging Face Transformers](https://github.com/huggingface/transformers)
- [TRL (SFTTrainer, DPO, etc.)](https://github.com/huggingface/trl)
- [Accelerate](https://github.com/huggingface/accelerate)
- [DeepSpeed](https://github.com/microsoft/DeepSpeed)
- [PyTorch FSDP docs](https://pytorch.org/docs/stable/fsdp.html)
- [Weights & Biases docs](https://docs.wandb.ai/)
- [MLflow docs](https://mlflow.org/docs/latest/index.html)
- [TensorBoard docs](https://www.tensorflow.org/tensorboard)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
