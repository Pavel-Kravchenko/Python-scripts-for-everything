---
name: ai-science-diffusion-generative-models
description: "Score matching, noise schedules, DDIM sampling, and DDRM inverse problems for diffusion generative models"
tool_type: python
primary_tool: NumPy
---

# Diffusion Generative Models

## Key Concepts

- **Forward (q)**: adds noise deterministically. **Reverse (p_θ)**: denoises via learned network
- **ᾱₜ vs αₜ**: αₜ = 1 − βₜ (single-step). ᾱₜ = ∏αₜ (cumulative). Forward closed form uses ᾱₜ
- **Cosine schedule** (Nichol & Dhariwal 2021): s=0.008 offset prevents ᾱ changing too fast near t=0
- **DDIM η**: η=0 → deterministic ODE sampler; η=1 → stochastic DDPM sampler

## Score Matching

Learn the score ∇_x log p(x) instead of p(x) directly. Denoising score matching:

$$\nabla_{x_t} \log p(x_t) \approx \frac{\hat{x}_0(x_t) - x_t}{\sigma_t^2}$$

**DDPM loss**: $\mathcal{L} = \mathbb{E}_{t, x_0, \varepsilon}[\|\varepsilon - \varepsilon_\theta(x_t, t)\|^2]$

## Noise Schedules

```python
import numpy as np

def linear_schedule(T=1000, beta_start=1e-4, beta_end=0.02):
    betas = np.linspace(beta_start, beta_end, T)
    alphas = 1 - betas
    return betas, alphas, np.cumprod(alphas)

def cosine_schedule(T=1000, s=0.008):
    t = np.arange(T + 1)
    f = np.cos((t / T + s) / (1 + s) * np.pi / 2) ** 2
    alpha_bars = f / f[0]
    betas = 1 - alpha_bars[1:] / alpha_bars[:-1]
    return np.clip(betas, 0, 0.999), 1 - np.clip(betas, 0, 0.999), alpha_bars[1:]
```

Forward diffusion: xₜ = √ᾱₜ · x₀ + √(1−ᾱₜ) · ε, where ε ~ N(0, I)

## DDIM Reverse Step

$$x_{t-1} = \sqrt{\bar\alpha_{t-1}} \hat{x}_0 + \sqrt{1 - \bar\alpha_{t-1} - \sigma_t^2} \cdot \varepsilon_\theta(x_t, t) + \sigma_t \varepsilon$$

Predicted clean image: $\hat{x}_0 = (x_t - \sqrt{1-\bar\alpha_t}\,\varepsilon_\theta) / \sqrt{\bar\alpha_t}$

```python
def ddim_step(xt, eps_pred, t, t_prev, alpha_bars, eta=0.0):
    ab_t = alpha_bars[t]
    ab_prev = alpha_bars[t_prev] if t_prev >= 0 else 1.0
    x0_pred = (xt - np.sqrt(1 - ab_t) * eps_pred) / (np.sqrt(ab_t) + 1e-9)
    x0_pred = np.clip(x0_pred, -1.5, 1.5)
    sigma = eta * np.sqrt((1 - ab_prev) / (1 - ab_t + 1e-9) * (1 - ab_t / ab_prev))
    noise = sigma * np.random.randn(*xt.shape) if eta > 0 else 0
    xt_prev = (np.sqrt(ab_prev) * x0_pred
               + np.sqrt(max(1 - ab_prev - sigma**2, 0)) * eps_pred
               + noise)
    return xt_prev, x0_pred
```

## DDRM: Inverse Problems via Diffusion

Degradation examples: denoising (H=I), inpainting (mask M), super-resolution (downsample D)

At each DDIM step, project x₀_pred onto measurement subspace via pseudoinverse of H (SVD):

$$\hat{x}_0^{\text{proj}} = \hat{x}_0 + V \cdot \frac{s^2}{s^2 + \sigma_y^2} \cdot \frac{U^T y - U^T H\hat{x}_0}{s}$$

```python
def ddrm_data_consistency(x0_pred, y_obs, H_type="identity", sigma_obs=0.3):
    if H_type == "identity":
        scale = 1 / (1 + sigma_obs**2)
        return x0_pred + scale * (y_obs - x0_pred)
    raise NotImplementedError(f"H_type '{H_type}' not implemented")
```

## Pitfalls

- Confusing forward (q, adds noise) with reverse (p_θ, denoises)
- Using αₜ where ᾱₜ is needed in the forward closed form
- Cosine schedule credit: Nichol & Dhariwal 2021, not the original Ho et al. 2020 DDPM

## Sources

- [DDRM (Kawar et al., 2022)](https://github.com/bahjat-kawar/ddrm)
