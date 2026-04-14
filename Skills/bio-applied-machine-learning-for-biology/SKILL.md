---
name: bio-applied-machine-learning-for-biology
description: Machine learning for bioinformatics — feature engineering for sequences, promoter classification, train/test splits, logistic regression, random forest, and bio-specific pitfalls
tool_type: python
primary_tool: Matplotlib
---

# Machine Learning for Biology

## ML Use Cases in Bioinformatics

| Application | Type | Key features | Tools |
|-------------|------|-------------|-------|
| Variant pathogenicity (CADD, REVEL) | Classification | Conservation, AF, protein impact | sklearn, XGBoost |
| Promoter / regulatory element | Classification | k-mer freq, CpG O/E, TATA box | sklearn |
| Cell type (scRNA-seq) | Clustering | Gene expression profiles | scanpy, sklearn |
| Drug response | Regression/Classification | Gene expression, mutations | sklearn |
| Protein function | Classification | Sequence features, GO terms | sklearn |
| Protein structure | Deep learning | Sequence + MSA | AlphaFold |

## Feature Engineering for Sequences

```python
from itertools import product

def kmer_frequencies(sequence, k=2):
    """Normalized k-mer frequencies. k=2 → 16 features, k=3 → 64 features."""
    sequence = sequence.upper()
    all_kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    counts = {km: 0 for km in all_kmers}
    total = 0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in counts:
            counts[kmer] += 1
            total += 1
    if total > 0:
        counts = {km: c/total for km, c in counts.items()}
    return counts

def extract_promoter_features(sequence, tss_pos=None):
    """Sequence-derived features for promoter classification."""
    seq = sequence.upper()
    n = len(seq)
    tss = tss_pos or n // 2
    features = {}

    features['gc_content'] = (seq.count('G') + seq.count('C')) / n

    for km, f in kmer_frequencies(seq, k=2).items():
        features[f'di_{km}'] = f

    n_c, n_g = seq.count('C'), seq.count('G')
    n_cpg = seq.count('CG')
    features['cpg_oe'] = (n_cpg * n) / (n_c * n_g) if n_c > 0 and n_g > 0 else 0

    upstream = seq[max(0, tss-50):tss]
    features['has_tata'] = 1 if 'TATAAA' in upstream else 0

    up = seq[max(0, tss-500):tss]
    dn = seq[tss:min(n, tss+500)]
    features['gc_upstream'] = (up.count('G')+up.count('C')) / max(len(up), 1)
    features['gc_downstream'] = (dn.count('G')+dn.count('C')) / max(len(dn), 1)

    return features
```

## Standard Workflow

```python
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score

# 1. Build feature matrix
feature_list = [extract_promoter_features(seq) for seq in sequences]
X = pd.DataFrame(feature_list).values
y = np.array(labels)

# 2. Stratified train/test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y)

# 3. Scale (required for logistic regression; not for tree models)
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_test_s = scaler.transform(X_test)

# 4a. Logistic regression (interpretable, linear boundary)
lr = LogisticRegression(max_iter=1000, random_state=42)
lr.fit(X_train_s, y_train)
print(classification_report(y_test, lr.predict(X_test_s)))

# 4b. Random forest (non-linear, feature importance, minimal tuning)
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)
print(classification_report(y_test, rf.predict(X_test)))

# 5. Cross-validation for reliable estimate
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(rf, X, y, cv=cv, scoring='roc_auc')
print(f"CV AUC: {scores.mean():.3f} ± {scores.std():.3f}")

# 6. Feature importance (random forest)
feat_imp = pd.Series(rf.feature_importances_, index=feature_cols).sort_values(ascending=False)
```

## Algorithm Decision Table

| Algorithm | Strengths | Weaknesses | Scale |
|-----------|-----------|------------|-------|
| Logistic Regression | Interpretable coefficients, fast | Linear boundary only | Fit first baseline |
| Random Forest | Non-linear, robust, feature importance | Less interpretable | Good default |
| SVM (RBF) | High-dim, small datasets | Slow to tune, no prob out of box | Try when RF overfits |
| Gradient Boosting (XGBoost) | Best tabular accuracy | More hyperparams | When RF not enough |
| Neural Net | Complex patterns, sequences | Data-hungry, black box | Large datasets |

## Pitfalls

- **Data leakage**: never fit scaler/imputer on the full dataset before splitting — fit only on training data, transform test separately
- **Class imbalance**: 99% negatives → 99% accuracy classifier that predicts all negative; use `stratify=y`, balanced class weights, or AUC-ROC instead of accuracy
- **Sequence similarity leakage**: train/test sequences from the same gene family or genomic region can share features — use chromosomal hold-out or CD-HIT clustering for proteins
- **k-mer curse of dimensionality**: k=6 → 4096 features for DNA; apply feature selection or PCA before linear models
- **Batch effects**: systematic differences between positive and negative set preparation confound results — balance data source, not just labels
- **Multiple testing**: apply Benjamini-Hochberg FDR when testing many features or models simultaneously
- **Overfitting small bio datasets**: typical promoter datasets have hundreds of examples; use cross-validation, not a single split
