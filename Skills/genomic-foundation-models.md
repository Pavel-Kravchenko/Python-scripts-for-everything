# genomic-foundation-models

DNA/RNA sequence foundation models: embeddings, fine-tuning, and regulatory prediction.

## Quick Reference

| Model | Parameters | Context | Best For |
|-------|-----------|---------|---------|
| DNABERT-2 | 117M | 512 bp | Classification tasks |
| Nucleotide Transformer | 2.5B | 6 kbp | General embeddings |
| HyenaDNA | 6.6M–1.6B | 1M bp | Long-range genomics |
| Evo | 7B | 131 kbp | Prokaryotic genomes, generation |
| Enformer | — | 196 kbp | Regulatory track prediction |

## Nucleotide Transformer Embeddings

```python
from transformers import AutoTokenizer, AutoModel
import torch
import numpy as np

model_name = 'InstaDeepAI/nucleotide-transformer-v2-500m-multi-species'
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name).eval()

sequences = ['ATCGATCGTAGCTAATCGATCG', 'GCTAGCTAGCTAGCTAGCTAG']

# Tokenize (6-mer tokens by default)
inputs = tokenizer(sequences, return_tensors='pt', padding=True, truncation=True, max_length=512)

with torch.no_grad():
    outputs = model(**inputs)

# Per-sequence embeddings (mean-pool over tokens)
embeddings = outputs.last_hidden_state.mean(dim=1).numpy()
# shape: (n_sequences, 1024)

# Per-residue embeddings for nucleotide-level tasks
token_embeddings = outputs.last_hidden_state.numpy()
```

## Fine-Tuning for Classification (Splice Sites)

```python
from transformers import AutoModelForSequenceClassification, TrainingArguments, Trainer
from datasets import Dataset
import torch

model = AutoModelForSequenceClassification.from_pretrained(
    'InstaDeepAI/nucleotide-transformer-v2-500m-multi-species',
    num_labels=2  # splice_site / non_splice_site
)

# Prepare dataset
train_dataset = Dataset.from_dict({
    'sequence': train_sequences,
    'label': train_labels
})

def tokenize(examples):
    return tokenizer(examples['sequence'], truncation=True, max_length=512, padding='max_length')

train_dataset = train_dataset.map(tokenize, batched=True)

# Training
args = TrainingArguments(
    output_dir='./nt_splice_finetuned',
    per_device_train_batch_size=8,
    num_train_epochs=3,
    learning_rate=1e-4,
    evaluation_strategy='epoch',
    fp16=True  # requires GPU
)

trainer = Trainer(model=model, args=args, train_dataset=train_dataset)
trainer.train()
```

## Enformer: Predict Regulatory Tracks

```python
# Using enformer-pytorch (community implementation)
# pip install enformer-pytorch

import torch
from enformer_pytorch import Enformer, from_pretrained

model = from_pretrained('EleutherAI/enformer-official-rough', target_length=-1)
model.eval()

# One-hot encode 196,608 bp sequence
# A=[1,0,0,0], C=[0,1,0,0], G=[0,0,1,0], T=[0,0,0,1], N=[0.25]*4
def one_hot_encode(seq):
    mapping = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1]}
    return torch.tensor([[mapping.get(c, [0.25]*4) for c in seq.upper()]], dtype=torch.float32)

sequence_196kbp = 'A' * 196608  # replace with real genomic sequence
one_hot = one_hot_encode(sequence_196kbp)  # shape: (1, 196608, 4)

with torch.no_grad():
    predictions = model(one_hot)
# predictions['human'] shape: (1, 896, 5313) → 5313 tracks at 128bp resolution
```

## In-Silico Mutagenesis (ISM) with Enformer

```python
import torch
import numpy as np

def ism_score(model, one_hot_seq, target_track=4799, center_pos=196608//2):
    """Score each possible SNP at center_pos for its effect on target track."""
    baseline = model(one_hot_seq)['human'][0, 448, target_track].item()
    scores = []
    for alt_nuc in range(4):  # A, C, G, T
        mut_seq = one_hot_seq.clone()
        mut_seq[0, center_pos, :] = 0
        mut_seq[0, center_pos, alt_nuc] = 1
        pred = model(mut_seq)['human'][0, 448, target_track].item()
        scores.append(pred - baseline)
    return scores  # [delta_A, delta_C, delta_G, delta_T]
```

## Common Pitfalls
- **Context length**: NT models use 6-mers → 512 tokens = ~3 kbp; for longer sequences use HyenaDNA
- **Tokenization**: DNABERT-2 uses BPE; NT uses overlapping 6-mers; different vocab sizes
- **Enformer input**: requires exactly 196,608 bp; pad or extract from genome
- **GPU memory**: NT-2.5B requires ~10GB VRAM; use half-precision (`model.half()`)
- **Species**: NT-multi trained on 850 species; NT-human better for human genomics

## Module
Tier 5 · Module 05 (Genomic Foundation Models)
