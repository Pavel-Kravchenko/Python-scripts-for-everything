---
name: bio-applied-machine-learning-for-biology
description: "Machine learning (ML) has become indispensable in modern bioinformatics: classifying protein functions, predicting drug responses, identifying regulatory elements, and analyzing single-cell data. This"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/07_Machine_Learning_for_Biology/01_machine_learning_for_biology.ipynb"
---

# Machine Learning for Biology

*Source: Course notebook `Tier_3_Applied_Bioinformatics/07_Machine_Learning_for_Biology/01_machine_learning_for_biology.ipynb`*

# Machine Learning for Biology

**Tier 3 -- Applied Bioinformatics**

Machine learning (ML) has become indispensable in modern bioinformatics: classifying protein functions, predicting drug responses, identifying regulatory elements, and analyzing single-cell data. This notebook introduces the core ML concepts and provides hands-on experience with scikit-learn on biological data.

**Prerequisites:** Tier 2 (Python, NumPy, pandas), statistics basics  
**Libraries:** `numpy`, `pandas`, `matplotlib`, `scikit-learn`

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

%matplotlib inline
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 12
np.random.seed(42)
```

---
## 1. Machine Learning in Bioinformatics: Use Cases

| Application | ML type | Input features | Output |
|-------------|---------|---------------|--------|
| Variant pathogenicity (CADD, REVEL) | Classification | Conservation, protein impact, allele frequency | Pathogenic / Benign |
| Gene expression prediction | Regression | Promoter sequence, histone marks | Expression level |
| Protein function prediction | Classification | Sequence features, domains, GO terms | Functional category |
| Drug response | Classification/Regression | Gene expression, mutations | Responder / IC50 |
| Cell type identification (scRNA-seq) | Clustering | Gene expression profiles | Cell type labels |
| Promoter recognition | Classification | k-mer frequencies, CpG content, TATA box | Promoter / Non-promoter |
| Protein structure (AlphaFold) | Deep learning | Sequence + MSA + templates | 3D coordinates |

---
## 2. Supervised Learning: Classification vs Regression

### Classification
The output is a **discrete category**: tumor vs. normal, protein family A vs. B vs. C.

### Regression
The output is a **continuous value**: gene expression level, drug IC50, protein stability.

Both follow the same workflow:
1. **Collect data** with known labels (training set)
2. **Engineer features** that capture relevant information
3. **Train a model** to map features to labels
4. **Evaluate** on held-out data
5. **Predict** on new, unlabeled data

---
## 3. Feature Engineering for Biological Data

Raw biological data (DNA sequences, protein structures) cannot be fed directly into most ML algorithms. We need to convert them into **numerical feature vectors**.

### 3.1 k-mer frequencies

For DNA/protein sequences, count the frequency of all subsequences of length k.

- **k=1**: nucleotide composition (4 features for DNA)
- **k=2**: dinucleotide frequencies (16 features for DNA)
- **k=3**: trinucleotide / codon frequencies (64 features for DNA)

```python
from itertools import product

def kmer_frequencies(sequence, k=2):
    """Compute normalized k-mer frequencies for a DNA sequence."""
    sequence = sequence.upper()
    # Generate all possible k-mers
    all_kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    counts = {km: 0 for km in all_kmers}
    
    total = 0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in counts:
            counts[kmer] += 1
            total += 1
    
    # Normalize to frequencies
    if total > 0:
        counts = {km: c / total for km, c in counts.items()}
    
    return counts

# Example
example_seq = 'ATCGATCGCGATCGATATCGCGCGATATATCGATCG'
freqs = kmer_frequencies(example_seq, k=2)
print("Dinucleotide frequencies:")
for km, f in sorted(freqs.items()):
    if f > 0:
        print(f"  {km}: {f:.3f}", end='  ')
print()
```

### 3.2 Physicochemical properties

For protein sequences, features can include:
- Amino acid composition (20 features)
- Molecular weight, isoelectric point
- Hydrophobicity profile
- Secondary structure propensities

### 3.3 Sequence-derived features for promoters

Drawing from the rice promoter dataset used in this course:
- GC content in different windows
- CpG observed/expected ratio
- TATA box presence/absence
- TFBS density and types
- Methylation levels
- SNP density

```python
def extract_promoter_features(sequence, tss_pos=None):
    """Extract a feature vector from a promoter sequence."""
    seq = sequence.upper()
    length = len(seq)
    if tss_pos is None:
        tss_pos = length // 2
    
    features = {}
    
    # Overall composition
    features['gc_content'] = (seq.count('G') + seq.count('C')) / length
    features['at_content'] = 1 - features['gc_content']
    
    # k-mer features (dinucleotides)
    di_freqs = kmer_frequencies(seq, k=2)
    for km, f in di_freqs.items():
        features[f'di_{km}'] = f
    
    # CpG observed/expected
    n_c = seq.count('C')
    n_g = seq.count('G')
    n_cpg = seq.count('CG')
    features['cpg_oe'] = (n_cpg * length) / (n_c * n_g) if (n_c > 0 and n_g > 0) else 0
    
    # TATA box presence
    upstream = seq[max(0, tss_pos-50):tss_pos]
    features['has_tata'] = 1 if 'TATAAA' in upstream else 0
    
    # GC content in upstream vs downstream windows
    up_500 = seq[max(0, tss_pos-500):tss_pos]
    down_500 = seq[tss_pos:min(length, tss_pos+500)]
    features['gc_upstream'] = (up_500.count('G') + up_500.count('C')) / max(len(up_500), 1)
    features['gc_downstream'] = (down_500.count('G') + down_500.count('C')) / max(len(down_500), 1)
    features['gc_ratio'] = features['gc_upstream'] / max(features['gc_downstream'], 0.01)
    
    return features

# Test on a random sequence
test_seq = ''.join(np.random.choice(list('ACGT'), 2000))
feats = extract_promoter_features(test_seq)
print(f"Number of features: {len(feats)}")
print("\nSample features:")
for k, v in list(feats.items())[:10]:
    print(f"  {k}: {v:.4f}")
```

---
## 4. Building a Biological Dataset: Promoter Classification

Let us create a realistic dataset for classifying DNA sequences as **promoter** or **non-promoter**. This is inspired by the rice promoter analysis project.

Promoter sequences will have:
- Higher GC content near the center
- Higher CpG density
- Occasional TATA boxes

Non-promoter sequences will have:
- Random genomic background composition

```python
def generate_promoter(length=500):
    """Generate a synthetic promoter sequence with biologically realistic features."""
    seq = []
    center = length // 2
    for i in range(length):
        dist = abs(i - center)
        gc_bias = 0.3 + 0.25 * np.exp(-0.5 * (dist / (length * 0.15)) ** 2)
        p = [0.5 * (1 - gc_bias), gc_bias / 2, gc_bias / 2, 0.5 * (1 - gc_bias)]
        seq.append(np.random.choice(list('ACGT'), p=p))
    
    # Add CpG enrichment
    for i in range(center - 100, center + 100):
        if i < length - 1 and np.random.random() < 0.15:
            seq[i] = 'C'
            seq[i + 1] = 'G'
    
    # Add TATA box with 30% probability
    if np.random.random() < 0.3:
        tata_pos = center - np.random.randint(25, 35)
        if tata_pos > 0 and tata_pos + 7 < length:
            for j, b in enumerate('TATAAAG'):
                seq[tata_pos + j] = b
    
    return ''.join(seq)

def generate_non_promoter(length=500):
    """Generate a random genomic (non-promoter) sequence."""
    weights = [0.29, 0.21, 0.21, 0.29]  # typical AT-rich vertebrate genome
    return ''.join(np.random.choice(list('ACGT'), size=length, p=weights))

# Generate dataset
np.random.seed(42)
n_pos = 300  # promoters
n_neg = 300  # non-promoters

sequences = []
labels = []

for _ in range(n_pos):
    sequences.append(generate_promoter())
    labels.append(1)

for _ in range(n_neg):
    sequences.append(generate_non_promoter())
    labels.append(0)

labels = np.array(labels)
print(f"Dataset: {len(sequences)} sequences")
print(f"  Promoters: {(labels == 1).sum()}")
print(f"  Non-promoters: {(labels == 0).sum()}")
```

```python
# Extract features for all sequences
feature_list = []
for seq in sequences:
    feature_list.append(extract_promoter_features(seq))

features_df = pd.DataFrame(feature_list)
features_df['label'] = labels

print(f"Feature matrix shape: {features_df.shape}")
print(f"\nFeature means by class:")
print(features_df.groupby('label')[['gc_content', 'cpg_oe', 'has_tata', 'di_CG']].mean().round(4))
```

```python
# Visualize feature distributions by class
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

plot_features = ['gc_content', 'cpg_oe', 'di_CG', 'gc_upstream']
plot_titles = ['GC Content', 'CpG Observed/Expected', 'CG Dinucleotide Freq', 'Upstream GC Content']

for ax, feat, title in zip(axes.flat, plot_features, plot_titles):
    prom_vals = features_df[features_df['label'] == 1][feat]
    non_vals = features_df[features_df['label'] == 0][feat]
    ax.hist(non_vals, bins=30, alpha=0.6, color='gray', label='Non-promoter', density=True)
    ax.hist(prom_vals, bins=30, alpha=0.6, color='steelblue', label='Promoter', density=True)
    ax.set_xlabel(title)
    ax.set_ylabel('Density')
    ax.legend()

plt.suptitle('Feature Distributions: Promoter vs Non-Promoter', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.show()
```

---
## 5. Train/Test Split and Cross-Validation

### Why we split data

A model must be evaluated on data it has **never seen during training**. Otherwise, we cannot distinguish genuine patterns from memorized noise (overfitting).

### Cross-validation

Instead of a single split, **k-fold cross-validation** rotates through k different train/test splits and averages the results. This gives a more reliable estimate of model performance.

```python
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold

# Prepare feature matrix and labels
feature_cols = [c for c in features_df.columns if c != 'label']
X = features_df[feature_cols].values
y = features_df['label'].values

# 80/20 train/test split (stratified to maintain class balance)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

print(f"Training set: {X_train.shape[0]} samples ({(y_train == 1).sum()} promoters, {(y_train == 0).sum()} non-promoters)")
print(f"Test set:     {X_test.shape[0]} samples ({(y_test == 1).sum()} promoters, {(y_test == 0).sum()} non-promoters)")
print(f"Features:     {X_train.shape[1]}")
```

---
## 6. Key Algorithms

### 6.1 Logistic Regression

Despite the name, logistic regression is a **classification** algorithm. It models the probability of class membership using the logistic function:

$$P(y=1|x) = \frac{1}{1 + e^{-(w^T x + b)}}$$

**Strengths:** Interpretable coefficients, fast, works well with many features.  
**Weaknesses:** Assumes linear decision boundary in feature space.

```python
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report

# Feature scaling (important for logistic regression)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train logistic regression
lr = LogisticRegression(max_iter=1000, random_state=42)
lr.fit(X_train_scaled, y_train)

y_pred_lr = lr.predict(X_test_scaled)
print("=== Logistic Regression ===")
print(f"Accuracy: {accuracy_score(y_test, y_pred_lr):.3f}")
print()
print(classification_report(y_test, y_pred_lr, target_names=['Non-promoter', 'Promoter']))
```

```python
# Feature importance from logistic regression coefficients
coef_df = pd.DataFrame({
    'feature': feature_cols,
    'coefficient': lr.coef_[0]
}).sort_values('coefficient', key=abs, ascending=False)

fig, ax = plt.subplots(figsize=(10, 6))
top_n = 15
top_feats = coef_df.head(top_n)
colors = ['steelblue' if c > 0 else 'coral' for c in top_feats['coefficient']]
ax.barh(range(top_n), top_feats['coefficient'], color=colors)
ax.set_yticks(range(top_n))
ax.set_yticklabels(top_feats['feature'])
ax.set_xlabel('Coefficient (positive = promoter-associated)')
ax.set_title('Top 15 Features -- Logistic Regression')
ax.invert_yaxis()
plt.tight_layout()
plt.show()
```

### 6.2 Random Forest

An **ensemble** of decision trees, each trained on a random subset of data and features. The final prediction is the majority vote.

**Strengths:** Handles non-linear relationships, robust to outliers, provides feature importance, minimal tuning needed.  
**Weaknesses:** Less interpretable than logistic regression, can overfit on very small datasets.

The source material for this course used Random Forest to classify wine varieties using chemical properties, achieving ~97% accuracy -- demonstrating the power of this approach.
