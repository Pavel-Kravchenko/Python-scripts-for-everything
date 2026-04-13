---
name: ai-science-diffusion-generative-models
description: "*Prerequisites: Module T5-01 (LLM Fine-tuning), Tier 3 Module 10 (Deep Learning), numpy.*"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb"
---

# Module T5-03: Diffusion & Generative Models

*Source: Course notebook `Tier_5_Modern_AI_for_Science/03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb`*

# Module T5-03: Diffusion & Generative Models

**Tier 5 — Modern AI for Science | Module 03**

*Prerequisites: Module T5-01 (LLM Fine-tuning), Tier 3 Module 10 (Deep Learning), numpy.*

GPU optional — all cells run on CPU with small examples.

---

**By the end you will be able to:**
- Explain score matching and the diffusion objective
- Implement linear and cosine noise schedules
- Run a DDIM sampling loop with an oracle denoiser
- Formulate image restoration as an inverse problem (DDRM)
- Implement SVD-based degradation operators
- Visualize score fields

**Attribution:** *Adapted from DDRM (Kawar et al., 2022), github.com/bahjat-kawar/ddrm. Uses public MNIST-style data.*

## How to work through this notebook

1. All cells run on CPU — this is intentional to keep diffusion mathematics accessible without requiring hardware.
2. Read the score-matching section (Section 1) before running code — the math connects forward and reverse processes.
3. Run the noise schedule plots before the DDIM section; understanding ᾱₜ is essential for the sampling loop.
4. The DDIM cells use an oracle denoiser (it knows x₀). Real networks learn this mapping. Focus on the loop logic.
5. DDRM (Sections 4–5) builds directly on DDIM; run those sections in order.

## Common sticking points

- **Forward vs reverse process direction**: forward (q) adds noise deterministically; reverse (p_θ) denoises stochastically using a learned neural network. Never confuse the two.
- **ᾱₜ vs αₜ**: αₜ = 1 − βₜ is the single-step retention. ᾱₜ = ∏αₜ is the cumulative retention (from t=0 to t). The forward diffusion closed form uses ᾱₜ, not αₜ.
- **Cosine schedule credit**: the cosine schedule was introduced by Nichol & Dhariwal (2021, "Improved DDPM"), not the original Ho et al. (2020) DDPM paper which used linear. The s=0.008 offset prevents ᾱ from changing too fast near t=0.
- **DDIM η parameter**: η=0 gives a deterministic ODE sampler (same latent → same output). η=1 recovers the stochastic DDPM sampler. η=0 allows latent interpolation; η=1 gives better diversity.
- **DDRM data consistency**: the projection step steers each denoising step toward observed measurements without refitting the network. It works because SVD decomposes H into independent singular-value channels.

## 1. Score-Based Generative Models

**Key insight:** Instead of learning p(x) directly, learn its *score* — the gradient of the log-density:

$$s_\theta(x) \approx \nabla_x \log p(x)$$

The score field points toward regions of higher probability. Starting from noise x_T ~ N(0, I), we can follow the score back to a data sample.

**Denoising Score Matching (DSM):** At noise level σₜ, a denoiser that predicts the clean image from the noisy one implicitly learns the score:

$$\nabla_{x_t} \log p(x_t) \approx \frac{\hat{x}_0(x_t) - x_t}{\sigma_t^2}$$

**DDPM objective:** Train ε_θ to predict the noise ε added to x₀:
$$\mathcal{L} = \mathbb{E}_{t, x_0, \varepsilon}\left[\|\varepsilon - \varepsilon_\theta(x_t, t)\|^2\right]$$

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({
    "figure.dpi": 120,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
print("Libraries loaded. All computations run on CPU.")
```

## 2. Noise Schedules

The noise schedule controls how quickly signal is destroyed as t increases from 0 (clean) to T (pure noise).

**Linear schedule (Ho et al., 2020 — original DDPM):** βₜ increases linearly from β₁ ≈ 10⁻⁴ to β_T ≈ 0.02. Works well for images but causes abrupt noise addition near t=0.

**Cosine schedule (Nichol & Dhariwal, 2021 — Improved DDPM):**

$$\bar\alpha_t = \frac{f(t)}{f(0)}, \quad f(t) = \cos^2\!\left(\frac{t/T + s}{1+s} \cdot \frac{\pi}{2}\right)$$

The small offset s=0.008 prevents ᾱ from changing too rapidly near t=0. The cosine schedule gives smoother signal destruction and generally improves sample quality.

**Critical quantity: ᾱₜ**
- ᾱₜ close to 1 → mostly signal (early timesteps)
- ᾱₜ close to 0 → mostly noise (late timesteps)
- Forward process closed form: xₜ = √ᾱₜ · x₀ + √(1−ᾱₜ) · ε, where ε ~ N(0, I)

```python
def linear_schedule(T=1000, beta_start=1e-4, beta_end=0.02):
    betas = np.linspace(beta_start, beta_end, T)
    alphas = 1 - betas
    alpha_bars = np.cumprod(alphas)
    return betas, alphas, alpha_bars

def cosine_schedule(T=1000, s=0.008):
    t = np.arange(T + 1)
    f = np.cos((t / T + s) / (1 + s) * np.pi / 2) ** 2
    alpha_bars = f / f[0]
    betas = 1 - alpha_bars[1:] / alpha_bars[:-1]
    betas = np.clip(betas, 0, 0.999)
    return betas, 1 - betas, alpha_bars[1:]

T = 1000
_, _, ab_linear  = linear_schedule(T)
_, _, ab_cosine  = cosine_schedule(T)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# ᾱₜ comparison
axes[0].plot(np.arange(T), ab_linear, label="Linear", color="steelblue")
axes[0].plot(np.arange(T), ab_cosine, label="Cosine", color="coral")
axes[0].set_xlabel("Timestep t"); axes[0].set_ylabel("ᾱₜ")
axes[0].set_title("Signal fraction ᾱₜ vs time"); axes[0].legend()

# Signal and noise amplitudes
axes[1].plot(np.arange(T), np.sqrt(ab_cosine), label="√ᾱₜ (signal)", color="steelblue")
axes[1].plot(np.arange(T), np.sqrt(1 - ab_cosine), label="√(1−ᾱₜ) (noise)", color="coral")
axes[1].set_xlabel("Timestep t"); axes[1].set_ylabel("Amplitude")
axes[1].set_title("Signal vs noise amplitudes (cosine)"); axes[1].legend()

# SNR
snr = ab_cosine / (1 - ab_cosine + 1e-9)
axes[2].semilogy(np.arange(T), snr, color="purple")
axes[2].set_xlabel("Timestep t"); axes[2].set_ylabel("SNR (log scale)")
axes[2].set_title("Signal-to-noise ratio (cosine schedule)")

plt.tight_layout(); plt.show()

# When does signal become < 1%?
t_threshold = np.argmax(ab_cosine < 0.01)
print(f"ᾱₜ < 0.01 at t = {t_threshold} (linear: {np.argmax(ab_linear < 0.01)})")
```

```python
# Forward diffusion: add noise to a simple 1D signal
n = 64
x0 = np.sin(np.linspace(0, 4*np.pi, n)) * 0.8  # clean "signal"

_, _, alpha_bars = cosine_schedule(T=1000)

def forward_diffuse(x0, t, alpha_bars, seed=0):
    rng = np.random.default_rng(seed)
    ab = alpha_bars[t]
    eps = rng.standard_normal(x0.shape)
    return np.sqrt(ab) * x0 + np.sqrt(1 - ab) * eps, eps

timesteps = [0, 100, 300, 500, 700, 999]
fig, axes = plt.subplots(1, len(timesteps), figsize=(15, 3))
for ax, t in zip(axes, timesteps):
    xt, _ = forward_diffuse(x0, t, alpha_bars)
    ax.plot(xt, lw=1.5, color="steelblue", alpha=0.8)
    ax.plot(x0, lw=1, color="red", alpha=0.4, ls="--", label="x₀" if t==0 else "")
    ax.set_title(f"t={t}\nᾱₜ={alpha_bars[t]:.3f}")
    ax.set_ylim(-3, 3); ax.set_xticks([]); ax.set_yticks([])

axes[0].legend(fontsize=8)
plt.suptitle("Forward diffusion: adding noise over time", y=1.05)
plt.tight_layout(); plt.show()
```

## 3. DDIM: Deterministic Sampling

DDPM sampling requires T=1000 steps (slow). DDIM reformulates as a non-Markovian process, enabling deterministic sampling in 20–50 steps.

**DDIM reverse step:**
$$x_{t-1} = \sqrt{\bar\alpha_{t-1}} \hat{x}_0 + \sqrt{1 - \bar\alpha_{t-1} - \sigma_t^2} \cdot \varepsilon_\theta(x_t, t) + \sigma_t \varepsilon$$

where η controls stochasticity: η=0 → fully deterministic (DDIM); η=1 → DDPM.

**Predicted clean image:**
$$\hat{x}_0 = \frac{x_t - \sqrt{1-\bar\alpha_t}\,\varepsilon_\theta(x_t, t)}{\sqrt{\bar\alpha_t}}$$

```python
def oracle_denoiser(xt, t, alpha_bars, x0_true):
    """
    Oracle denoiser: predicts the true noise ε from xt.
    In practice, this is replaced by a neural network.
    """
    ab = alpha_bars[t]
    # Tweedie's formula: E[x0 | xt] = (xt - sqrt(1-ab)*eps) / sqrt(ab)
    # We know x0_true, so we compute the implied eps
    eps_true = (xt - np.sqrt(ab) * x0_true) / np.sqrt(1 - ab + 1e-9)
    return eps_true

def ddim_step(xt, eps_pred, t, t_prev, alpha_bars, eta=0.0):
    """One DDIM reverse step."""
    ab_t    = alpha_bars[t]
    ab_prev = alpha_bars[t_prev] if t_prev >= 0 else 1.0

    # Predicted x0
    x0_pred = (xt - np.sqrt(1 - ab_t) * eps_pred) / (np.sqrt(ab_t) + 1e-9)
    x0_pred = np.clip(x0_pred, -1.5, 1.5)

    # DDIM step
    sigma = eta * np.sqrt((1 - ab_prev) / (1 - ab_t + 1e-9) * (1 - ab_t / ab_prev))
    noise = sigma * np.random.randn(*xt.shape) if eta > 0 else 0
    xt_prev = (np.sqrt(ab_prev) * x0_pred
               + np.sqrt(max(1 - ab_prev - sigma**2, 0)) * eps_pred
               + noise)
    return xt_prev, x0_pred

# Demo: DDIM reverse on our 1D signal
_, _, alpha_bars = cosine_schedule(T=1000)
x0_true = np.sin(np.linspace(0, 4*np.pi, 64)) * 0.8

# Start from noisy image at T=999
xT, _ = forward_diffuse(x0_true, 999, alpha_bars, seed=1)

# DDIM sampling with T_sample steps
T_sample = 20
timesteps = np.linspace(999, 1, T_sample).astype(int)

xt = xT.copy()
trajectories = [xT.copy()]
for i, t in enumerate(timesteps[:-1]):
    t_prev = timesteps[i+1]
    eps_pred = oracle_denoiser(xt, t, alpha_bars, x0_true)
    xt, x0_pred = ddim_step(xt, eps_pred, t, t_prev, alpha_bars, eta=0.0)
    trajectories.append(xt.copy())

fig, axes = plt.subplots(1, 5, figsize=(15, 3))
steps_to_show = [0, 5, 10, 15, T_sample-1]
for ax, step_idx in zip(axes, steps_to_show):
    ax.plot(trajectories[step_idx], color="steelblue", lw=1.5)
    ax.plot(x0_true, color="red", lw=1, ls="--", alpha=0.5, label="x₀ true")
    t_label = timesteps[step_idx] if step_idx < len(timesteps) else 0
    ax.set_title(f"Step {step_idx}/{T_sample}\nt={t_label}")
    ax.set_ylim(-3, 3); ax.set_xticks([]); ax.set_yticks([])

axes[0].legend(fontsize=8)
plt.suptitle("DDIM reverse: recovering signal from noise (20 steps, oracle denoiser)", y=1.05)
plt.tight_layout(); plt.show()

mse = ((trajectories[-1] - x0_true)**2).mean()
print(f"Final MSE (oracle denoiser): {mse:.6f}")
```

## 4. Inverse Problems (DDRM)

Many scientific imaging tasks are *inverse problems*: given a degraded measurement y, recover the original signal x₀.

**Examples:**
| Degradation | H matrix | y = Hx₀ + noise |
|---|---|---|
| Denoising | Identity H=I | y = x₀ + σn |
| Inpainting | Mask M | y = M⊙x₀ |
| Super-resolution | Downsampling D | y = Dx₀ |
| Colorization | Luminance L | y = L(x₀) |

**DDRM approach:** At each DDIM step, project x₀_pred onto the measurement subspace using the pseudoinverse of H (via SVD), then continue sampling.

$$\hat{x}_0^{\text{proj}} = \hat{x}_0 + V \cdot \frac{s^2}{s^2 + \sigma_y^2} \cdot \frac{U^T y - U^T H\hat{x}_0}{s}$$

where H = UΣVᵀ is the SVD of H.

```python
# DDRM: denoising as inverse problem
# H = I (identity), y = x0 + sigma_obs * noise
sigma_obs = 0.3
rng = np.random.default_rng(5)
y_obs = x0_true + sigma_obs * rng.standard_normal(x0_true.shape)

def ddrm_data_consistency(x0_pred, y_obs, H_type="identity", sigma_obs=0.3):
    """Apply DDRM data-consistency step for H=I (denoising)."""
    if H_type == "identity":
        # For H=I, SVD is trivial: U=I, s=1, V=I
        # Projection: x0_proj = x0_pred + (s²/(s²+σ²)) * (y - x0_pred)
        scale = 1 / (1 + sigma_obs**2)
        return x0_pred + scale * (y_obs - x0_pred)
    else:
        raise NotImplementedError(f"H_type '{H_type}' not implemented")

# DDIM + DDRM sampling
xt = xT.copy()
traj_ddrm = [xT.copy()]

for i, t in enumerate(timesteps[:-1]):
    t_prev = timesteps[i+1]
    eps_pred = oracle_denoiser(xt, t, alpha_bars, x0_true)
    xt_candidate, x0_pred = ddim_step(xt, eps_pred, t, t_prev, alpha_bars, eta=0.0)
    # Data consistency
    x0_proj = ddrm_data_consistency(x0_pred, y_obs, sigma_obs=sigma_obs)
    # Recompute noise from projected x0
    ab_t = alpha_bars[t_prev]
    eps_proj = (xt_candidate - np.sqrt(ab_t) * x0_proj) / np.sqrt(1 - ab_t + 1e-9)
    xt = xt_candidate
    traj_ddrm.append(xt.copy())

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(x0_true,      color="black", lw=2, label="x₀ (true)", ls="-")
ax.plot(y_obs,        color="red", alpha=0.5, lw=1, label=f"y (noisy, σ={sigma_obs})", ls="--")
ax.plot(traj_ddrm[-1], color="steelblue", lw=2, label="DDRM reconstruction", alpha=0.9)
ax.plot(trajectories[-1], color="green", lw=1.5, ls=":", label="DDIM (no data consistency)")
ax.set_xlabel("Signal position"); ax.set_ylabel("Amplitude")
ax.set_title("DDRM: denoising via diffusion + data consistency")
ax.legend(frameon=False); plt.tight_layout(); plt.show()

mse_naive = ((y_obs - x0_true)**2).mean()
mse_ddrm  = ((traj_ddrm[-1] - x0_true)**2).mean()
print(f"MSE (noisy observation): {mse_naive:.4f}")
print(f"MSE (DDRM):              {mse_ddrm:.4f}")
print(f"DDRM improvement: {(mse_naive - mse_ddrm)/mse_naive:.1%}")
```
