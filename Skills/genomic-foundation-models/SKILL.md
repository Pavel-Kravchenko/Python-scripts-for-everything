---
name: genomic-foundation-models
description: "DNA/RNA sequence foundation models: embeddings, fine-tuning, and regulatory prediction."
tool_type: python
primary_tool: HuggingFace Transformers
---

# genomic-foundation-models

## Model Landscape

| Model | Context | Best For |
|-------|---------|---------|
| DNABERT-2 | ~3 kbp (BPE tokens) | Short-window classification (promoter/enhancer/splice) |
| Nucleotide Transformer | kb-scale (6-mer tokens) | Transfer learning across species |
| HyenaDNA | up to 1M bp | Distal regulatory context |
| Evo | ~100 kb+ | Prokaryotic sequence design (autoregressive) |
| Enformer | 196,608 bp | Multi-track regulatory signal (5313 tracks, 128 bp bins) |
| Borzoi | 524 kb | RNA-seq coverage at 32 bp bins |
| SpliceAI | local window | Clinical splice variant delta scores |
| AlphaGenome | up to 1M bp | Unified variant effect (expression + splicing + chromatin + contacts) |

**Routing:** General embeddings: DNABERT-2 / NT / HyenaDNA. Regulatory prediction: Enformer / Borzoi / AlphaGenome. Splice scoring: SpliceAI. Sequence design: Evo.

### Complementary Protein Structure Models

| Model | Best For |
|------|---------|
| AlphaFold 2 | Protein monomer structure |
| AlphaFold 3 | Complexes (proteins + nucleic acids + ligands) |
| RoseTTAFold | Open academic alternative |

## Nucleotide Transformer Embeddings

```python
from transformers import AutoTokenizer, AutoModel
import torch

model_name = 'InstaDeepAI/nucleotide-transformer-v2-500m-multi-species'
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name).eval()

inputs = tokenizer(sequences, return_tensors='pt', padding=True, truncation=True, max_length=512)
with torch.no_grad():
    outputs = model(**inputs)

# Per-sequence: mean-pool over tokens -> (n_sequences, 1024)
embeddings = outputs.last_hidden_state.mean(dim=1).numpy()
# Per-residue: outputs.last_hidden_state -> (n_sequences, seq_len_tokens, 1024)
```

## Fine-Tuning for Classification

```python
from transformers import AutoModelForSequenceClassification, TrainingArguments, Trainer

model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)
args = TrainingArguments(
    output_dir='./finetuned', per_device_train_batch_size=8,
    num_train_epochs=3, learning_rate=1e-4, evaluation_strategy='epoch', fp16=True
)
trainer = Trainer(model=model, args=args, train_dataset=tokenized_dataset)
trainer.train()
```

## Enformer: Regulatory Track Prediction

```python
from enformer_pytorch import from_pretrained
import torch

model = from_pretrained('EleutherAI/enformer-official-rough', target_length=-1).eval()

def one_hot_encode(seq):
    mapping = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1]}
    return torch.tensor([[mapping.get(c, [0.25]*4) for c in seq.upper()]], dtype=torch.float32)

# Input: exactly 196,608 bp one-hot -> output: (1, 896, 5313) tracks at 128 bp resolution
with torch.no_grad():
    predictions = model(one_hot_encode(sequence_196kbp))
```

## In-Silico Mutagenesis (ISM)

```python
def ism_score(model, one_hot_seq, target_track=4799, center_pos=196608//2):
    baseline = model(one_hot_seq)['human'][0, 448, target_track].item()
    scores = []
    for alt_nuc in range(4):
        mut_seq = one_hot_seq.clone()
        mut_seq[0, center_pos, :] = 0
        mut_seq[0, center_pos, alt_nuc] = 1
        pred = model(mut_seq)['human'][0, 448, target_track].item()
        scores.append(pred - baseline)
    return scores  # [delta_A, delta_C, delta_G, delta_T]
```

## Pitfalls
- **Context length**: NT 6-mers with 512 tokens is only ~3 kbp; use HyenaDNA for longer
- **Tokenization mismatch**: DNABERT-2 uses BPE; NT uses overlapping 6-mers -- different vocab
- **Enformer input**: requires exactly 196,608 bp; pad or extract from genome
- **GPU memory**: NT-2.5B needs ~10 GB VRAM; use `model.half()`
- **Species fit**: NT-multi (850 species) vs NT-human -- pick the right checkpoint
- **Scope**: SpliceAI is splice-only; Enformer/Borzoi/AlphaGenome are broader regulatory models
- **Domain confusion**: AlphaFold/RoseTTAFold predict protein structure, not genomic tracks
