---
name: machine-learning-bio
description: Classical machine learning for biological tabular and sequence-derived features.
---

## When to Use

Use this atomic skill for focused work on **machine-learning-bio** without bundling unrelated topics.

## Quick Reference

This skill was split from `ml-deep-learning-bio.md` to keep topics independent and self-contained.

## Core Patterns

Use the parent material below as the source reference, then keep implementations specific to this topic.

## Source Reference (from merged skill)

---
name: ml-deep-learning-bio
description: Machine learning and deep learning for biological data — scikit-learn classifiers, CNNs for genomic sequences, protein language models
---

# ML & Deep Learning for Biology

## When to Use
- Classifying sequences (promoter vs. non-promoter, pathogenic vs. benign variant)
- Predicting protein function, secondary structure, or expression levels
- Dimensionality reduction / clustering of scRNA-seq or expression profiles
- Motif detection in DNA/protein sequences with CNNs
- Transfer learning with ESM-2, DNABERT-2, or AlphaFold embeddings

---

## Quick Reference

### Algorithm selection

| Situation | Use |
|-----------|-----|
| Tabular features, interpretability required | Logistic Regression |
| Non-linear, robust baseline, feature importance | Random Forest |
| High-dimensional features (k-mers, expression) | SVM (RBF kernel) |
| Small dataset, no tuning | k-NN |
| Raw sequences, motif detection | 1D CNN |
| Long-range dependencies, full-length sequences | Bi-LSTM or Transformer |
| Per-residue/per-position prediction | Bi-LSTM or Transformer encoder |
| Denoising, latent space, generation | VAE |
| Small data + relevant pretrained model available | Transfer learning (ESM-2, DNABERT-2) |

### Classical ML vs. deep learning

- **Classical ML** (RF, SVM, LR): tabular features, small datasets (<1000 samples), interpretability critical
- **Deep learning**: raw sequences/images, large datasets (>10k samples), hierarchical patterns

### Metrics summary

| Metric | Use when |
|--------|----------|
| Accuracy | Balanced classes |
| F1 | Imbalanced classes (e.g. pathogenic variants ~5% positive) |
| ROC-AUC | Ranking quality, threshold-independent |
| Precision-Recall AUC | Highly imbalanced classes |

---

## Key Patterns

### Feature engineering for sequences

**k-mer frequencies** (DNA/protein):
- k=1: nucleotide composition (4 features)
- k=2: dinucleotide frequencies (16 features)
- k=3: trinucleotide / codon frequencies (64 features)

**Promoter-specific features**: GC content, CpG observed/expected ratio, TATA box presence, regional GC asymmetry

**Protein features**: amino acid composition (20), molecular weight, isoelectric point, hydrophobicity profile

**One-hot encoding** for CNN input: shape `(N, 4, L)` for DNA — channels-first for `Conv1d`

### Regularization for biological data ("small n, large p")
- **Lasso (L1)**: use when few features drive prediction (sparse model, e.g. key marker genes)
- **Ridge (L2)**: use when many features contribute small effects
- **ElasticNet**: default when unsure

---

## Code Templates

### Scikit-learn full pipeline

```python
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import f1_score

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

pipelines = {
    'LR':  Pipeline([('sc', StandardScaler()), ('clf', LogisticRegression(max_iter=1000))]),
    'RF':  Pipeline([('clf', RandomForestClassifier(n_estimators=200, max_depth=10, n_jobs=-1))]),
    'SVM': Pipeline([('sc', StandardScaler()), ('clf', SVC(kernel='rbf', probability=True))]),
    'kNN': Pipeline([('sc', StandardScaler()), ('clf', KNeighborsClassifier(n_neighbors=7))]),
}

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
for name, pipe in pipelines.items():
    scores = cross_val_score(pipe, X, y, cv=cv, scoring='roc_auc')
    print(f"{name:5s}  AUC {scores.mean():.3f} ± {scores.std():.3f}")

# Hyperparameter search
param_grid = {'clf__n_estimators': [100, 200, 500], 'clf__max_depth': [5, 10, None]}
gs = GridSearchCV(pipelines['RF'], param_grid, cv=5, scoring='f1', n_jobs=-1)
gs.fit(X_train, y_train)
print(gs.best_params_, f1_score(y_test, gs.predict(X_test)))
```

### k-mer feature extraction

```python
from itertools import product

def kmer_frequencies(sequence, k=2):
    sequence = sequence.upper()
    all_kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    counts = {km: 0 for km in all_kmers}
    total = 0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in counts:
            counts[kmer] += 1
            total += 1
    return {km: c / total for km, c in counts.items()} if total > 0 else counts
```

### One-hot encoding for CNN input

```python
import numpy as np, torch

def one_hot_encode(sequences, alphabet='ACGT'):
    """Returns tensor (N, len(alphabet), L) — Conv1d channels-first format."""
    mapping = {c: i for i, c in enumerate(alphabet)}
    n, L = len(sequences), len(sequences[0])
    enc = np.zeros((n, len(alphabet), L), dtype=np.float32)
    for i, seq in enumerate(sequences):
        for j, c in enumerate(seq):
            if c in mapping:
                enc[i, mapping[c], j] = 1.0
    return torch.FloatTensor(enc)
```

### PyTorch training loop (standard)

```python
import torch, torch.nn as nn, torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
train_dl = DataLoader(TensorDataset(X_train_t, y_train_t), batch_size=32, shuffle=True)

model = MyModel().to(device)
optimizer = optim.Adam(model.parameters(), lr=1e-3)
criterion = nn.BCELoss()

for epoch in range(100):
    model.train()
    for xb, yb in train_dl:
        xb, yb = xb.to(device), yb.to(device)
        optimizer.zero_grad()
        loss = criterion(model(xb), yb)
        loss.backward()
        optimizer.step()

model.eval()
with torch.no_grad():
    preds = model(X_test_t.to(device)).cpu()
    acc = ((preds > 0.5).float() == y_test_t).float().mean()
```

### 1D CNN for DNA sequences

```python
class SequenceCNN(nn.Module):
    def __init__(self, seq_length=200):
        super().__init__()
        self.conv_layers = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=8, padding=3), nn.ReLU(), nn.MaxPool1d(4),
            nn.Conv1d(32, 64, kernel_size=6, padding=2), nn.ReLU(), nn.MaxPool1d(4),
        )
        self.classifier = nn.Sequential(
            nn.Linear(64, 64), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(64, 1), nn.Sigmoid()
        )
    def forward(self, x):          # x: (batch, 4, L)
        x = self.conv_layers(x)    # -> (batch, 64, L')
        x = x.mean(dim=2)          # global average pooling -> (batch, 64)
        return self.classifier(x)
```

### Bidirectional LSTM for sequences

```python
class SequenceLSTM(nn.Module):
    def __init__(self, input_dim=4, hidden_dim=64, num_layers=2, dropout=0.3):
        super().__init__()
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers,
                            batch_first=True, bidirectional=True, dropout=dropout)
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim * 2, 64), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 1), nn.Sigmoid()
        )
    def forward(self, x):                           # x: (batch, 4, L)
        x = x.transpose(1, 2)                       # -> (batch, L, 4)
        _, (h_n, _) = self.lstm(x)
        combined = torch.cat([h_n[-2], h_n[-1]], dim=1)  # fwd + bwd last hidden
        return self.classifier(combined)
```

### Transformer encoder for sequences

```python
class TransformerClassifier(nn.Module):
    def __init__(self, input_dim=4, d_model=64, nhead=4, num_layers=2, max_len=200):
        super().__init__()
        self.input_proj = nn.Linear(input_dim, d_model)
        self.pos_embedding = nn.Parameter(torch.randn(1, max_len, d_model) * 0.02)
        enc_layer = nn.TransformerEncoderLayer(
            d_model, nhead, dim_feedforward=128, dropout=0.1, batch_first=True)
        self.transformer = nn.TransformerEncoder(enc_layer, num_layers)
        self.classifier = nn.Sequential(
            nn.Linear(d_model, 32), nn.ReLU(), nn.Linear(32, 1), nn.Sigmoid())
    def forward(self, x):
        x = x.transpose(1, 2)                               # (batch, L, 4)
        x = self.input_proj(x) + self.pos_embedding[:, :x.size(1)]
        x = self.transformer(x).mean(dim=1)                 # global avg pool
        return self.classifier(x)
```

### VAE for gene expression / latent space

```python
class VAE(nn.Module):
    def __init__(self, input_dim, hidden_dim=128, latent_dim=10):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim), nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim // 2), nn.ReLU())
        self.fc_mu     = nn.Linear(hidden_dim // 2, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim // 2, latent_dim)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2), nn.ReLU(),
            nn.Linear(hidden_dim // 2, hidden_dim), nn.ReLU(),
            nn.Linear(hidden_dim, input_dim))
    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h)
    def reparameterize(self, mu, logvar):
        return mu + torch.exp(0.5 * logvar) * torch.randn_like(logvar)
    def forward(self, x):
        mu, logvar = self.encode(x)
        return self.decoder(self.reparameterize(mu, logvar)), mu, logvar

def vae_loss(recon, x, mu, logvar, beta=1.0):
    recon_loss = nn.functional.mse_loss(recon, x, reduction='sum')
    kl_loss    = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return recon_loss + beta * kl_loss
```

### DNA data augmentation

```python
def augment_dna_sequence(seq):
    complement = str.maketrans('ACGT', 'TGCA')
    rc = seq[::-1].translate(complement)          # reverse complement
    mutated = list(seq)
    for _ in range(3):                            # 3 random point mutations
        mutated[np.random.randint(len(seq))] = np.random.choice(list('ACGT'))
    return [rc, ''.join(mutated)]
```

### Saliency map (CNN interpretability)

```python
def compute_saliency(model, x, device):
    """Which sequence positions most influence the prediction."""
    x = x.unsqueeze(0).requires_grad_(True)
    model.eval()
    model(x.to(device)).sum().backward()
    return x.grad.data.abs().squeeze(0).sum(dim=0).cpu().numpy()  # shape: (L,)
```

---

## Common Pitfalls

- **Data leakage**: fit scaler on train only — `.fit_transform(X_train)`, `.transform(X_test)` — never the reverse
- **Evaluating on training data**: cross-validation must wrap the full pipeline including the scaler
- **Class imbalance**: use `stratify=y` in split; `class_weight='balanced'` in LR/SVM; report F1/AUC, not accuracy
- **RF does not need scaling**: trees are scale-invariant; scaling required for LR, SVM, kNN
- **CNN input shape**: PyTorch `Conv1d` expects `(batch, channels, length)` — one-hot DNA is `(N, 4, L)`, not `(N, L, 4)`
- **Bi-LSTM input**: needs `(batch, seq_len, features)` — transpose from one-hot with `.transpose(1, 2)`
- **Overfitting in DL**: biology has small n, large p — always use dropout + early stopping; watch validation loss
- **ROC-AUC = 0.5**: model is random — check for data leakage or label errors before tuning
- **k-NN with many features**: degrades badly with >50 features; apply PCA first

---

## Pre-trained Models for Transfer Learning

| Model | Domain | Embedding dim | Access |
|-------|--------|--------------|--------|
| ESM-2 (Meta AI) | Proteins | 320–5120 | `pip install fair-esm` |
| ProtT5 (Rostlab) | Proteins | 1024 | HuggingFace `Rostlab/prot_t5_xl_uniref50` |
| DNABERT-2 | DNA (multi-species) | 768 | HuggingFace `zhihan1996/DNABERT-2-117M` |
| Nucleotide Transformer | DNA | 512–1024 | HuggingFace `InstaDeepAI/nucleotide-transformer-*` |

Use embeddings as frozen features for small datasets; fine-tune last layers when >1000 labeled examples are available.

---

## Related Skills
- `numpy-pandas-wrangling` — expression matrix preprocessing before ML
- `biostatistics-r` — statistical testing to complement ML feature selection
- `biopython-databases` — sequence retrieval and feature extraction pipelines


## Related Skills

- `machine-learning-bio` (this file)
- `ml-deep-learning-bio` (legacy merged skill)
