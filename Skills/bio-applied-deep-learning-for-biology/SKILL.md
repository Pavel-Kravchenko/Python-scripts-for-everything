---
name: bio-applied-deep-learning-for-biology
description: PyTorch deep learning for biological sequences — when to use DL vs classical ML, CNN architecture for motif detection, one-hot encoding, and the standard training loop
tool_type: python
primary_tool: NumPy
---

# Deep Learning for Biology

## Classical ML vs Deep Learning Decision Table

| Criterion | Classical ML | Deep Learning |
|-----------|-------------|---------------|
| Sample size | 100s–1000s | 10,000+ (or transfer learning) |
| Feature engineering | Manual (k-mers, physicochemical) | Learned automatically |
| Input type | Tabular features | Raw sequences / images |
| Interpretability | High (feature importance) | Lower (requires SHAP/attention) |
| Training time | Minutes | Hours–days |
| Hardware | CPU | GPU recommended |

Use classical ML (RF, SVM) when: tabular features, small dataset, interpretability required.
Use DL when: raw sequence/image input, large dataset, hierarchical patterns, pre-trained models available.

## PyTorch Core Patterns

```python
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Autograd
x = torch.tensor([1.0, 2.0, 3.0], requires_grad=True)
loss = (x ** 2).sum()
loss.backward()       # x.grad == 2*x

# NumPy <-> Tensor (shares memory for CPU tensors)
t = torch.from_numpy(np_array.astype(np.float32))
arr = tensor.detach().numpy()
```

## Model Definition

```python
class SequenceClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)
```

## Standard Training Loop

```python
model = SequenceClassifier(input_dim=16).to(device)
criterion = nn.BCELoss()           # binary; use CrossEntropyLoss for multiclass
optimizer = optim.Adam(model.parameters(), lr=1e-3)
train_loader = DataLoader(TensorDataset(X_train_t, y_train_t), batch_size=32, shuffle=True)

for epoch in range(100):
    model.train()
    for X_batch, y_batch in train_loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)
        optimizer.zero_grad()             # 1. reset gradients
        preds = model(X_batch)            # 2. forward pass
        loss = criterion(preds, y_batch)  # 3. compute loss
        loss.backward()                   # 4. backprop
        optimizer.step()                  # 5. update weights

model.eval()
with torch.no_grad():
    test_preds = model(X_test_t.to(device)).cpu()
    acc = ((test_preds > 0.5).float() == y_test_t).float().mean()
```

## One-Hot Encoding for DNA Sequences

```python
def one_hot_encode(sequences, alphabet='ACGT'):
    """Returns tensor shape (N, 4, L) — Conv1d expects (batch, channels, length)."""
    mapping = {c: i for i, c in enumerate(alphabet)}
    n, L = len(sequences), len(sequences[0])
    enc = np.zeros((n, len(alphabet), L), dtype=np.float32)
    for i, seq in enumerate(sequences):
        for j, c in enumerate(seq):
            if c in mapping:
                enc[i, mapping[c], j] = 1.0
    return torch.FloatTensor(enc)
```

## 1D CNN for Motif Detection

```python
class SeqCNN(nn.Module):
    """Conv1D architecture: detect motifs -> pool -> classify."""
    def __init__(self, n_filters=32, kernel_size=8):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(4, n_filters, kernel_size),   # 4 channels = ACGT
            nn.ReLU(),
            nn.MaxPool1d(4),                         # position invariance
            nn.Conv1d(n_filters, n_filters * 2, kernel_size // 2),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1)                  # global pooling
        )
        self.fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(n_filters * 2, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.fc(self.conv(x))
```

## Dinucleotide Features (Classical Input)

```python
def dinucleotide_features(sequences):
    """16 dinucleotide frequencies per sequence — fast baseline."""
    dinucs = [a+b for a in 'ACGT' for b in 'ACGT']
    idx = {d: i for i, d in enumerate(dinucs)}
    X = np.zeros((len(sequences), 16))
    for i, seq in enumerate(sequences):
        for j in range(len(seq) - 1):
            d = seq[j:j+2]
            if d in idx:
                X[i, idx[d]] += 1
        X[i] /= (len(seq) - 1)
    return X
```

## Pitfalls

- **`model.eval()` + `torch.no_grad()`**: both required for inference — `eval()` disables dropout/batchnorm; `no_grad()` prevents gradient tracking (memory/speed)
- **`optimizer.zero_grad()` before backward**: forgetting this accumulates gradients across batches and produces wrong updates
- **Conv1d input shape**: expects `(batch, channels, length)`; one-hot DNA is `(N, 4, L)` — transpose from `(N, L, 4)` if needed
- **Vanishing gradients**: avoid sigmoid/tanh in hidden layers; use ReLU; add BatchNorm for deep networks
- **Class imbalance**: use `pos_weight` in `BCELoss` or weighted sampler — biological datasets are often highly imbalanced (e.g., 1% binding sites vs background)
- **Data leakage**: split train/test BEFORE any scaling/normalization; fit scaler on train only, apply to test
- **GPU memory**: move model AND data to same device; check `tensor.device` when debugging cryptic errors
