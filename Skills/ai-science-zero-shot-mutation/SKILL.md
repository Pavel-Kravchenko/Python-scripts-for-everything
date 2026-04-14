---
name: ai-science-zero-shot-mutation
description: Zero-Shot Mutation Effect Prediction with NumPy
tool_type: python
primary_tool: NumPy
---

# Zero-Shot Mutation Effect Prediction

## Key Distinctions

- **Masked vs autoregressive**: ESM-1v uses MLM — mask position i, compare log-probs. Autoregressive models (ProtGPT2) compare full sequence log-likelihoods.
- **Single-site vs epistatic**: ESM-1v computes single-site marginals independently. For epistatic questions use ESM-IF1 or DMS-fine-tuned models.
- **ΔΔG vs ΔlogP**: ESM scores track evolutionary fitness, not thermodynamic stability — correlated but not identical.
- **Benchmark first**: check Spearman ρ on ProteinGym DMS assays for your protein family before trusting scores.

```python
import numpy as np
import pandas as pd

np.random.seed(13)
AA = list('ACDEFGHIKLMNPQRSTVWY')
```

## Core Formula

For a masked LM like ESM-1v, the zero-shot score for wt → mutant at position i:

$$\Delta_{i, \text{wt}\to\text{mut}} = \log P_\theta(\text{mut} \mid x_{\setminus i}) - \log P_\theta(\text{wt} \mid x_{\setminus i})$$

where $x_{\setminus i}$ is the sequence with position $i$ masked. Negative values mean the mutant is less likely under the model's evolutionary prior — a signal for deleterious effect.

ESM-1v trained on ~250M sequences: masking position i lets the model reconstruct it from residue covariance. Evolutionarily tolerated substitutions get high probability; disruptive ones get low.

## Real ESM-1v Usage (requires GPU + `fair-esm`)

```python
import esm
import torch
import torch.nn.functional as F

# Load ESM-1v (5 ensemble members available use model 1 as default)
model, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
batch_converter = alphabet.get_batch_converter()
model.eval()

def score_mutations(sequence: str, positions: list[int]) -> dict:
    \"\"\"
    Compute ESM-1v zero-shot ΔlogP for all mutations at given positions.
    positions: 1-indexed residue positions to score
    \"\"\"
    data = [("seq", sequence)]
    _, _, tokens = batch_converter(data)

    scores = {}
    with torch.no_grad():
        for pos in positions:  # pos is 1-indexed
            # Mask position i (tokens are 1-indexed due to <cls>)
            masked_tokens = tokens.clone()
            masked_tokens[0, pos] = alphabet.mask_idx  # pos+1 for <cls> offset

            logits = model(masked_tokens, repr_layers=[])["logits"]  # (1, L, vocab)
            log_probs = F.log_softmax(logits[0, pos], dim=-1)

            wt_aa = sequence[pos - 1]
            wt_idx = alphabet.get_idx(wt_aa)

            pos_scores = {}
            for mut_aa in "ACDEFGHIKLMNPQRSTVWY":
                mut_idx = alphabet.get_idx(mut_aa)
                delta = (log_probs[mut_idx] - log_probs[wt_idx]).item()
                pos_scores[f"{wt_aa}{pos}{mut_aa}"] = delta
            scores.update(pos_scores)
    return scores
```

## Ensemble Averaging

ESM-1v has 5 independently trained checkpoints. Averaging across all 5 significantly improves Spearman ρ with experimental fitness.

```python
# Load all 5 ESM-1v models
models_and_alphabets = [
    esm.pretrained.esm1v_t33_650M_UR90S_1(),
    esm.pretrained.esm1v_t33_650M_UR90S_2(),
    esm.pretrained.esm1v_t33_650M_UR90S_3(),
    esm.pretrained.esm1v_t33_650M_UR90S_4(),
    esm.pretrained.esm1v_t33_650M_UR90S_5(),
]
# Average the delta scores from all 5 models for each variant
```

## Toy Proxy (no GPU needed)

```python
HYDRO = set('AILMFWVY')
POLAR = set('STNQ')
CHARGED = set('KRHDE')

def aa_group(a):
    if a in HYDRO: return 'hydro'
    if a in POLAR: return 'polar'
    if a in CHARGED: return 'charged'
    return 'other'

def toy_zero_shot_delta(wt: str, mut: str) -> float:
    # penalize physicochemical group switches
    return 0.25 if aa_group(wt) == aa_group(mut) else -0.55
```

## Score All Single Substitutions

```python
wt_seq = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ'

records = []
for i, wt in enumerate(wt_seq):
    for mut in AA:
        if mut == wt:
            continue
        delta = toy_zero_shot_delta(wt, mut)
        records.append({'pos': i + 1, 'wt': wt, 'mut': mut, 'delta_ll': delta})

df = pd.DataFrame(records)
df.sort_values('delta_ll').head(8)
```

## Benchmark Against Pseudo-Fitness

```python
noise = np.random.normal(0, 0.08, size=len(df))
df['fitness'] = 0.6 + 0.5 * df['delta_ll'] + noise

corr = df[['delta_ll', 'fitness']].corr(method='spearman').iloc[0, 1]
print('Spearman correlation (delta vs fitness):', round(float(corr), 3))
```

## Clinical Triage Table

Combine zero-shot effect with rarity and structural confidence proxies.

```python
subset = df.sample(12, random_state=1).copy()
subset['rarity'] = np.random.uniform(0.2, 1.0, size=len(subset))
subset['plddt_proxy'] = np.random.uniform(55, 95, size=len(subset))

subset['priority'] = (
    0.5 * (-subset['delta_ll']) +
    0.3 * subset['rarity'] +
    0.2 * (subset['plddt_proxy'] / 100.0)
)

subset.sort_values('priority', ascending=False)[['pos', 'wt', 'mut', 'delta_ll', 'rarity', 'plddt_proxy', 'priority']]
```

## References

- [ESM repository (ESM-1v/ESM-2)](https://github.com/facebookresearch/esm)
- [ProteinGym benchmark](https://proteingym.org/)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
