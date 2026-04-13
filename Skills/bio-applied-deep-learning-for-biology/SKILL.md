---
name: bio-applied-deep-learning-for-biology
description: "Deep learning has driven breakthroughs across biology: AlphaFold solved the protein structure prediction problem, transformer-based protein language models capture evolutionary information from raw se"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/10_Deep_Learning_for_Biology/01_deep_learning_for_biology.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, pytorch 2.2+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Deep Learning for Biology

*Source: Course notebook `Tier_3_Applied_Bioinformatics/10_Deep_Learning_for_Biology/01_deep_learning_for_biology.ipynb`*


**Tier 3 -- Applied Bioinformatics**

Deep learning has driven breakthroughs across biology: AlphaFold solved the protein structure prediction problem, transformer-based protein language models capture evolutionary information from raw sequences, and variational autoencoders reveal latent structure in single-cell data. This notebook builds from classical ML foundations (covered in Module 07) and introduces the deep learning architectures most relevant to bioinformatics.

**Prerequisites:** Module 07 (Machine Learning for Biology), NumPy, pandas, matplotlib  
**Libraries:** `torch` (PyTorch), `numpy`, `pandas`, `matplotlib`

```python
# Installation (run once)
# For CPU-only PyTorch (sufficient for this notebook):
#   pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
# For GPU-enabled PyTorch (requires CUDA):
#   pip install torch torchvision
# On Google Colab, PyTorch is pre-installed.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

%matplotlib inline
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 12
np.random.seed(42)
torch.manual_seed(42)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")
print(f"PyTorch version: {torch.__version__}")
```python

---
## 1. From Classical ML to Deep Learning

### When Do You Need Deep Learning?

Classical ML (random forests, SVMs, logistic regression) remains the right choice for many biological problems, especially when:
- You have **tabular features** (gene expression, k-mer counts, physicochemical properties)
- Datasets are **small** (hundreds to low thousands of samples)
- **Interpretability** is critical

Deep learning excels when:
- Input is **raw sequence/structure/image** data (no hand-crafted features needed)
- Datasets are **large** (tens of thousands or more samples)
- The problem has **hierarchical patterns** (motifs within motifs)
- You can leverage **pre-trained models** via transfer learning

| Criterion | Classical ML | Deep Learning |
|-----------|-------------|---------------|
| Sample size needed | 100s--1000s | 1000s--millions |
| Feature engineering | Manual (domain knowledge) | Learned automatically |
| Interpretability | High (feature importance) | Lower (requires special methods) |
| Raw sequences/images | Needs preprocessing | Native input |
| Training time | Minutes | Hours--days |
| Hardware | CPU | GPU recommended |

### Neural Network Fundamentals

A neural network is a series of **layers** that transform input data through learnable linear transformations followed by nonlinear **activation functions**.

**The Perceptron (single neuron):**

$$y = \sigma(\mathbf{w} \cdot \mathbf{x} + b)$$

where $\mathbf{w}$ are weights, $b$ is a bias, and $\sigma$ is an activation function.

**Common activation functions:**
- **Sigmoid:** $\sigma(z) = \frac{1}{1 + e^{-z}}$ -- outputs in (0, 1), used for binary classification output
- **ReLU:** $\sigma(z) = \max(0, z)$ -- default choice for hidden layers, avoids vanishing gradients
- **Tanh:** $\sigma(z) = \tanh(z)$ -- outputs in (-1, 1), centered around zero

**Backpropagation** computes the gradient of the loss with respect to every weight using the chain rule, allowing gradient descent to update the weights.

**Universal Approximation Theorem (intuition):** A neural network with a single hidden layer containing enough neurons can approximate any continuous function to arbitrary precision. In practice, *deeper* networks learn hierarchical representations more efficiently than very wide shallow ones.

```python
# Visualize activation functions
z = np.linspace(-5, 5, 200)

sigmoid = 1 / (1 + np.exp(-z))
relu = np.maximum(0, z)
tanh = np.tanh(z)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for ax, vals, name in zip(axes, [sigmoid, relu, tanh], ['Sigmoid', 'ReLU', 'Tanh']):
    ax.plot(z, vals, linewidth=2)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5)
    ax.set_title(name)
    ax.set_xlabel('z')
    ax.set_ylabel('activation(z)')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```python

---
## 2. Building Neural Networks with PyTorch

### Tensors and Autograd

PyTorch tensors are like NumPy arrays but can track gradients for automatic differentiation.

```python
# Tensors: creating and basic operations
x = torch.tensor([1.0, 2.0, 3.0], requires_grad=True)
y = x ** 2 + 2 * x + 1  # elementwise quadratic
loss = y.sum()

# Backpropagation: compute dy/dx
loss.backward()

print(f"x     = {x.data}")
print(f"y     = {y.data}")
print(f"dy/dx = {x.grad}")  # dy/dx = 2x + 2
print(f"Expected: {2 * x.data + 2}")
```python

```python
# Converting between NumPy and PyTorch
np_array = np.array([[1, 2], [3, 4]], dtype=np.float32)
tensor = torch.from_numpy(np_array)
back_to_np = tensor.numpy()

print(f"NumPy shape: {np_array.shape}, dtype: {np_array.dtype}")
print(f"Tensor shape: {tensor.shape}, dtype: {tensor.dtype}")
print(f"Shares memory: {np.shares_memory(np_array, back_to_np)}")
```python

### Defining Models with nn.Module

Every PyTorch model inherits from `nn.Module`. You define the layers in `__init__` and the forward pass in `forward`. Backpropagation through the layers is handled automatically.

```python
class SimpleClassifier(nn.Module):
    """A feedforward network for binary classification."""

    def __init__(self, input_dim, hidden_dim=64):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.network(x)

# Inspect the model
model = SimpleClassifier(input_dim=16)
print(model)
total_params = sum(p.numel() for p in model.parameters())
print(f"\nTotal parameters: {total_params:,}")
```python

### Case Study: Promoter vs Non-Promoter Classification

We build a binary classifier using **dinucleotide frequencies** as features. This parallels the classical ML approach from Module 07 but uses a neural network.

```python
# Generate synthetic promoter / non-promoter sequences
def generate_sequences(n_per_class=500, seq_length=200):
    sequences, labels = [], []
    for _ in range(n_per_class):
        # Promoter: enriched in CG dinucleotides, TATA-like motifs
        seq = []
        for j in range(seq_length):
            if np.random.random() < 0.15:
                seq.append(np.random.choice(['C', 'G']))
            elif 80 <= j <= 90 and np.random.random() < 0.6:
                seq.extend(list('TATAAA'[:min(6, seq_length - j)]))
            else:
                seq.append(np.random.choice(['A', 'T', 'G', 'C'], p=[0.2, 0.2, 0.35, 0.25]))
        sequences.append(''.join(seq[:seq_length]))
        labels.append(1)

        # Non-promoter: uniform composition
        seq = ''.join(np.random.choice(['A', 'T', 'G', 'C'], size=seq_length))
        sequences.append(seq)
        labels.append(0)
    return sequences, np.array(labels)


def dinucleotide_features(sequences):
    """Extract 16 dinucleotide frequencies from each sequence."""
    nucleotides = ['A', 'C', 'G', 'T']
    dinucs = [a + b for a in nucleotides for b in nucleotides]
    features = np.zeros((len(sequences), 16))
    for i, seq in enumerate(sequences):
        total = len(seq) - 1
        for j in range(total):
            dinuc = seq[j:j+2]
            if dinuc in dinucs:
                features[i, dinucs.index(dinuc)] += 1
        features[i] /= total  # normalize to frequencies
    return features, dinucs


sequences, labels = generate_sequences(500)
X, feature_names = dinucleotide_features(sequences)
print(f"Dataset: {X.shape[0]} sequences, {X.shape[1]} features")
print(f"Labels: {np.sum(labels == 1)} promoters, {np.sum(labels == 0)} non-promoters")
```python

```python
# The standard PyTorch training loop
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Split and scale
X_train, X_test, y_train, y_test = train_test_split(
    X, labels, test_size=0.2, random_state=42, stratify=labels
)
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Convert to tensors
X_train_t = torch.FloatTensor(X_train)
y_train_t = torch.FloatTensor(y_train).unsqueeze(1)
X_test_t = torch.FloatTensor(X_test)
y_test_t = torch.FloatTensor(y_test).unsqueeze(1)

# DataLoader for mini-batch training
train_dataset = TensorDataset(X_train_t, y_train_t)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)

# Model, loss, optimizer
model = SimpleClassifier(input_dim=16, hidden_dim=64).to(device)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Training loop
train_losses = []
for epoch in range(100):
    model.train()
    epoch_loss = 0
    for X_batch, y_batch in train_loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)

        optimizer.zero_grad()           # 1. Reset gradients
        predictions = model(X_batch)    # 2. Forward pass
        loss = criterion(predictions, y_batch)  # 3. Compute loss
        loss.backward()                 # 4. Backpropagation
        optimizer.step()                # 5. Update weights

        epoch_loss += loss.item() * X_batch.size(0)

    train_losses.append(epoch_loss / len(train_dataset))

# Evaluate
model.eval()
with torch.no_grad():
    test_preds = model(X_test_t.to(device)).cpu()
    test_acc = ((test_preds > 0.5).float() == y_test_t).float().mean()

print(f"Final training loss: {train_losses[-1]:.4f}")
print(f"Test accuracy: {test_acc:.4f}")
```python

```python
# Plot training loss
plt.figure(figsize=(8, 4))
plt.plot(train_losses, linewidth=2)
plt.xlabel('Epoch')
plt.ylabel('Binary Cross-Entropy Loss')
plt.title('Training Loss -- Promoter Classifier')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```python

---
## 3. Convolutional Neural Networks for Sequences

### One-Hot Encoding of Biological Sequences

Instead of extracting hand-crafted features (k-mer frequencies), we can feed **raw sequences** into a CNN. Each position becomes a 4-dimensional one-hot vector (for DNA) or a 20-dimensional vector (for proteins).

```python
Sequence: A  C  G  T  A
A:        1  0  0  0  1
C:        0  1  0  0  0
G:        0  0  1  0  0
T:        0  0  0  1  0
```python

A 1D convolution slides a **filter** (kernel) across the sequence. Filters learn to detect motifs -- just like k-mers, but the network discovers the relevant patterns itself.

```python
def one_hot_encode(sequences, alphabet='ACGT'):
    """One-hot encode a list of sequences.

    Returns tensor of shape (N, C, L) where C = len(alphabet), L = seq_length.
    This matches PyTorch Conv1d expected input: (batch, channels, length).
    """
    mapping = {c: i for i, c in enumerate(alphabet)}
    n = len(sequences)
    seq_len = len(sequences[0])
    encoded = np.zeros((n, len(alphabet), seq_len), dtype=np.float32)
    for i, seq in enumerate(sequences):
        for j, char in enumerate(seq[:seq_len]):
            if char in mapping:
                encoded[i, mapping[char], j] = 1.0
    return torch.FloatTensor(encoded)


# Encode our sequences
X_encoded = one_hot_encode(sequences)
print(f"Encoded shape: {X_encoded.shape}  (samples, channels=4, length=200)")

# Visualize one-hot encoding for a short segment
fig, ax = plt.subplots(figsize=(12, 2))
segment = X_encoded[0, :, :30].numpy()
ax.imshow(segment, aspect='auto', cmap='Blues', interpolation='nearest')
ax.set_yticks(range(4))
ax.set_yticklabels(['A', 'C', 'G', 'T'])
ax.set_xlabel('Position')
ax.set_title('One-hot encoding (first 30 bp of a promoter)')
plt.tight_layout()
plt.show()
```python

### 1D CNN Architecture for TF Binding Site Prediction

The architecture follows a proven pattern for biological sequence classification:

**Conv1D** (detect motifs) --> **ReLU** (nonlinearity) --> **MaxPool** (position invariance) --> **Dense** (classify)

Multiple convolutional layers can capture motif combinations and spacing patterns.

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
