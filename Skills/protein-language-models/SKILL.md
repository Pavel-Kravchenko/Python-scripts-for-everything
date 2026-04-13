---
name: protein-language-models
description: "ESM2 embeddings, ESMFold structure prediction, zero-shot mutation scoring, and protein design."
tool_type: python
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pytorch 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# protein-language-models

ESM2 embeddings, ESMFold structure prediction, zero-shot mutation scoring, and protein design.

## Quick Reference

| Model | Parameters | Best For |
|-------|-----------|---------|
| ESM2-8M | 8M | Fast embeddings, resource-constrained |
| ESM2-650M | 650M | Good quality, fits T4 GPU |
| ESM2-3B | 3B | High quality, requires A100 |
| ESM2-15B | 15B | State-of-the-art, multi-GPU |
| ESMFold | 690M | Fast structure prediction (no MSA) |
| ESM-1v | 650M×5 | Zero-shot mutation scoring |
| ESM-IF1 | 142M | Inverse folding / sequence design |

## ESM2 Sequence Embeddings

```python
import esm
import torch
import numpy as np

# Load model and alphabet
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model.eval()
batch_converter = alphabet.get_batch_converter()

# Prepare sequences
data = [
    ('protein_A', 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD'),
    ('protein_B', 'MKTLLLTLVVVTIVCLDLGYTPETRFLLKFNKAVIVAGTNTK'),
]
batch_labels, batch_strs, batch_tokens = batch_converter(data)

# Extract embeddings (layer 33)
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

# Per-residue embeddings: (batch, seq_len, 1280)
token_rep = results['representations'][33]

# Per-sequence embeddings: mean-pool over sequence length
seq_embeddings = token_rep[:, 1:-1, :].mean(dim=1).numpy()  # exclude BOS/EOS tokens
print(seq_embeddings.shape)  # (2, 1280)
```

## ESMFold Structure Prediction

```python
import esm
import torch

# Load ESMFold (requires ~16GB RAM, GPU recommended)
model = esm.pretrained.esmfold_v1()
model = model.eval().to('cuda')

sequence = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK'

with torch.no_grad():
    output = model.infer_pdb(sequence)

# Save PDB
with open('esmfold_structure.pdb', 'w') as f:
    f.write(output)

# Access pLDDT confidence scores
with torch.no_grad():
    output_dict = model.infer(sequence)
plddt = output_dict['plddt'].squeeze().cpu().numpy()  # per-residue confidence
print(f'Mean pLDDT: {plddt.mean():.2f}')  # > 70 = confident, > 90 = very confident
```

## Zero-Shot Mutation Scoring (ESM-1v)

```python
import esm
import torch

# Load ESM-1v (ensemble of 5 models for better accuracy)
model, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
model.eval()
batch_converter = alphabet.get_batch_converter()

def score_mutation(wt_seq, position, wt_aa, mut_aa):
    """
    Score effect of mutation using masked marginal log-likelihood.
    position: 0-indexed
    Returns: log P(mut | context) - log P(wt | context)
    """
    data = [('protein', wt_seq)]
    _, _, batch_tokens = batch_converter(data)
    
    # Mask the target position
    masked_tokens = batch_tokens.clone()
    masked_tokens[0, position + 1] = alphabet.mask_idx  # +1 for BOS token
    
    with torch.no_grad():
        logits = model(masked_tokens)['logits']
    
    log_probs = torch.nn.functional.log_softmax(logits[0, position + 1], dim=-1)
    
    wt_idx = alphabet.get_idx(wt_aa)
    mut_idx = alphabet.get_idx(mut_aa)
    
    return (log_probs[mut_idx] - log_probs[wt_idx]).item()

# Example: score effect of A2G mutation in BRCA1
wt_sequence = 'MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNACKESMLCTHSSLNFFPVSLNLNPFQNNRNQLQNELREQLKLRQLEMDLNRFLSEYRSSMSLNHLENSSAAQLKLMQQKELNQIFQELNFQNQNQNQNQNQ'
delta_llr = score_mutation(wt_sequence, 1, 'D', 'E')  # D2E mutation
print(f'Delta log-likelihood: {delta_llr:.3f}')  # negative = destabilizing
```

## Contact Map Prediction from Attention

```python
import esm
import torch
import matplotlib.pyplot as plt

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model.eval()
batch_converter = alphabet.get_batch_converter()

data = [('protein', 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD')]
_, _, tokens = batch_converter(data)

with torch.no_grad():
    results = model(tokens, repr_layers=[33], return_contacts=True)

contacts = results['contacts'][0]  # (L, L) contact map
plt.matshow(contacts.numpy(), cmap='RdYlBu_r')
plt.title('ESM2 Contact Prediction')
plt.colorbar()
```

## ESM-IF1 Inverse Folding

```python
import esm
import torch
import biotite.structure.io as bsio

model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
model.eval()

# Load PDB structure
structure = esm.inverse_folding.util.load_structure('protein.pdb', 'A')
coords, seq = esm.inverse_folding.util.extract_coords_from_structure(structure)

# Sample diverse sequences for the given backbone
sampled_seqs = model.sample(coords, temperature=1.0, partial_seq=None)
# temperature: 1.0=diverse, 0.1=conservative (close to native)
```

## Common Pitfalls
- **Padding**: ESM pads to max_len in batch; use attention mask for accurate mean-pooling
- **Token offset**: ESM adds `<cls>` at position 0; residue i is at token index i+1
- **ESMFold speed**: ~0.5s per sequence on GPU vs hours for AlphaFold2; accuracy slightly lower
- **Zero-shot scoring**: works best for single mutations; epistasis not captured
- **Sequence length**: ESM2 max = 1024 tokens; for longer proteins, use windowing

## Module
Tier 5 · Module 06 (Protein Language Models)
