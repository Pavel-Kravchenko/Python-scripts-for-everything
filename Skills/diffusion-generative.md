---
name: diffusion-generative
description: Score-based generative models, DDIM sampling, linear and cosine noise schedules, inverse imaging problems (DDRM), SVD degradation operators, scientific image restoration
---

## When to Use

Use this skill when:
- Implementing or understanding DDIM/DDPM sampling
- Computing linear or cosine noise schedules
- Formulating image restoration as an inverse problem (DDRM)
- Visualizing score fields (∇ₓ log p(x))
- Understanding diffusion for scientific imaging (cryo-EM, MRI)

## Quick Reference

| Symbol | Meaning | Formula |
|---|---|---|
| βₜ | Noise added at step t | schedule: linear or cosine |
| αₜ | Signal retained at step t | αₜ = 1 − βₜ |
| ᾱₜ | Cumulative signal fraction | ᾱₜ = ∏_{s=1}^{t} αₛ |
| xₜ | Noisy image at step t | xₜ = √ᾱₜ x₀ + √(1−ᾱₜ) ε, ε∼N(0,I) |
| ε_θ | Denoiser network prediction | predicts noise ε from xₜ |
| x₀_pred | Predicted clean image | x₀_pred = (xₜ − √(1−ᾱₜ) ε_θ) / √ᾱₜ |
| DDIM step | Deterministic reverse | xₜ₋₁ = √ᾱₜ₋₁ x₀_pred + √(1−ᾱₜ₋₁) ε_θ |

## Key Patterns

**Pattern 1: Linear noise schedule**
```python
import numpy as np

def linear_schedule(T=1000, beta_start=1e-4, beta_end=0.02):
    betas = np.linspace(beta_start, beta_end, T)
    alphas = 1 - betas
    alpha_bars = np.cumprod(alphas)
    return betas, alphas, alpha_bars
```

**Pattern 2: Cosine schedule (improved)**
```python
def cosine_schedule(T=1000, s=0.008):
    t = np.arange(T + 1)
    f = np.cos((t / T + s) / (1 + s) * np.pi / 2) ** 2
    alpha_bars = f / f[0]
    betas = 1 - alpha_bars[1:] / alpha_bars[:-1]
    return np.clip(betas, 0, 0.999), 1 - betas, alpha_bars[1:]
```

**Pattern 3: Forward diffusion (add noise)**
```python
def forward_diffuse(x0, t, alpha_bars):
    """Add noise to x0 at timestep t."""
    ab = alpha_bars[t]
    eps = np.random.randn(*x0.shape)
    xt = np.sqrt(ab) * x0 + np.sqrt(1 - ab) * eps
    return xt, eps
```

**Pattern 4: DDIM reverse step**
```python
def ddim_step(xt, eps_pred, t, t_prev, alpha_bars, eta=0.0):
    """One DDIM reverse step. eta=0 → deterministic."""
    ab_t    = alpha_bars[t]
    ab_prev = alpha_bars[t_prev] if t_prev >= 0 else 1.0
    x0_pred = (xt - np.sqrt(1 - ab_t) * eps_pred) / np.sqrt(ab_t)
    x0_pred = np.clip(x0_pred, -1, 1)
    sigma   = eta * np.sqrt((1 - ab_prev) / (1 - ab_t) * (1 - ab_t / ab_prev))
    noise   = sigma * np.random.randn(*xt.shape) if eta > 0 else 0
    xt_prev = np.sqrt(ab_prev) * x0_pred + np.sqrt(1 - ab_prev - sigma**2) * eps_pred + noise
    return xt_prev, x0_pred
```

**Pattern 5: SVD degradation (DDRM)**
```python
def apply_degradation(x, H):
    """Apply linear degradation H to image x. H: (m, n) matrix."""
    return H @ x.flatten()

def pseudoinverse_step(y_obs, x0_pred, H, sigma_obs=0.05):
    """DDRM data-consistency projection."""
    # Singular value decomposition of H
    U, s, Vt = np.linalg.svd(H, full_matrices=False)
    # Project onto measurement subspace
    Hy  = U.T @ y_obs
    Hx0 = U.T @ (H @ x0_pred.flatten())
    # Weighted update
    scale = s**2 / (s**2 + sigma_obs**2)
    x0_consistent = Vt.T @ (Vt @ x0_pred.flatten() + scale * (Hy - Hx0) / s)
    return x0_consistent.reshape(x0_pred.shape)
```

## Code Templates

**Template 1: DDIM sampling loop**
```python
def ddim_sample(oracle_denoiser, shape, T=100, T_sample=20, schedule="cosine"):
    if schedule == "cosine":
        betas, alphas, alpha_bars = cosine_schedule(T)
    else:
        betas, alphas, alpha_bars = linear_schedule(T)

    timesteps = np.linspace(T-1, 0, T_sample).astype(int)
    xt = np.random.randn(*shape)

    for i, t in enumerate(timesteps[:-1]):
        t_prev = timesteps[i+1]
        eps_pred = oracle_denoiser(xt, t)
        xt, x0_pred = ddim_step(xt, eps_pred, t, t_prev, alpha_bars)

    return xt
```

## Common Pitfalls

- **Tensor shape:** always use explicit shapes; `betas` must be (T,), not (T, 1) — broadcasting errors are silent
- **alpha_bar indexing:** alpha_bars[T-1] → most noisy; alpha_bars[0] → least noisy; check your indexing convention
- **Clipping x0_pred:** without clipping to [-1, 1], artifacts accumulate over DDIM steps
- **SVD conditioning:** near-zero singular values in H make pseudoinverse unstable; use truncated SVD with threshold
- **Cosine schedule for images:** cosine schedule gives smoother noise addition than linear; prefer it for image data

## Related Skills

- `ml-deep-learning-bio` — PyTorch and neural network training
- `llm-finetuning` — separate generative paradigm (autoregressive vs diffusion)
