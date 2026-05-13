# Bifurcation Forecasting for Flutter: Algorithm Description

*Based on Lim & Epureanu (2011), Ghadami & Epureanu (2016, 2017), Ghadami, Cesnik & Epureanu (2018), Riso, Ghadami, Cesnik & Epureanu (2020), Riso, Cesnik & Epureanu (2021, 2022).*

---

## 1. Overview and Physical Basis

Bifurcation forecasting predicts flutter onset and post-flutter limit-cycle oscillation (LCO) amplitudes using only pre-flutter free-decay time histories. No system model is required for the primary method. The method exploits **critical slowing down (CSD)**: as a control parameter $\boldsymbol{\mu}$ (e.g.\ flow speed, stiffness) approaches the flutter boundary, the system recovers from perturbations more slowly. Quantifying this slow-down allows the bifurcation diagram to be reconstructed by extrapolation rather than direct simulation.

**Key quantities:**
- $\boldsymbol{\mu} \in \mathbb{R}^{N_p}$: vector of $N_p$ control parameters.
- $r$: scalar response amplitude (distance from equilibrium, measured at fixed oscillation phase).
- $\lambda(\boldsymbol{\mu}, r) = \dot{r}/r = \mathrm{d}(\ln r)/\mathrm{d}t$: **recovery rate** (units: s$^{-1}$). Negative in pre-flutter regime; zero on the bifurcation diagram.

**Governing principle.** The amplitude evolution near a Hopf bifurcation can be written

$$\dot{r} = r\!\left[\alpha_0(r) + \boldsymbol{\alpha}_1(r)^T(\boldsymbol{\mu} - \tilde{\boldsymbol{\mu}}_c) + (\boldsymbol{\mu} - \tilde{\boldsymbol{\mu}}_c)^T \boldsymbol{\alpha}_2(r)(\boldsymbol{\mu} - \tilde{\boldsymbol{\mu}}_c) + \mathrm{HOT}\right]$$

where $\tilde{\boldsymbol{\mu}}_c$ is a point on the flutter boundary and the coefficients $\alpha_0, \boldsymbol{\alpha}_1, \boldsymbol{\alpha}_2$ are polynomial functions of $r$ that are independent of $\boldsymbol{\mu}$. This Taylor expansion is in parameter space only; the dynamics is **not** linearised in amplitude.

The recovery rate is therefore

$$\lambda(\boldsymbol{\mu}, r) \approx \alpha_0(r) + \boldsymbol{\alpha}_1(r)^T(\boldsymbol{\mu} - \tilde{\boldsymbol{\mu}}_c) + (\boldsymbol{\mu} - \tilde{\boldsymbol{\mu}}_c)^T \boldsymbol{\alpha}_2(r)(\boldsymbol{\mu} - \tilde{\boldsymbol{\mu}}_c)$$

Points on the bifurcation diagram satisfy $\lambda(\boldsymbol{\mu}, r) = 0$.

---

## 2. Method Assumptions

The following assumptions must hold for the method to be valid:

1. The system is smooth and varies smoothly with $\boldsymbol{\mu}$ near the bifurcation.
2. The system exhibits CSD as $\boldsymbol{\mu}$ approaches the flutter boundary.
3. Exactly one mode (the **bifurcating mode**) approaches instability; its associated eigenvalue pair crosses the imaginary axis. The centre space is two-dimensional.
4. The bifurcating mode is excited in the transient responses.
5. Each transient response converges to the original equilibrium state (perturbations must not push the system to a different attractor).
6. Control-parameter samples $\{\boldsymbol{\mu}_l\}$ must be close enough to the flutter boundary that measurable CSD is present.

---

## 3. Data Requirements

Collect $N_s$ free-decay time histories in the pre-flutter regime:

- Choose $N_s$ control-parameter vectors $\{\boldsymbol{\mu}_l\}_{l=1}^{N_s}$, all strictly below the flutter boundary, sampled in a neighbourhood of the expected flutter point.
  - Use **Latin Hypercube Sampling (LHS)** within a small parameter box around the expected flutter region.
  - For a **first-order** (linear) Taylor expansion with $N_p$ parameters: minimum $N_s = 1 + N_p$ transients.
  - For a **second-order** (quadratic) Taylor expansion: minimum $N_s = 1 + N_p + N_p(N_p+1)/2$.
  - More transients improve robustness to noise; overdetermined systems are solved by least squares.
- At each $\boldsymbol{\mu}_l$, apply a **large-amplitude perturbation** (e.g.\ a 1-cosine gust, an impulse, or a non-zero initial condition) and record the **free-decay response** $y(t; \boldsymbol{\mu}_l)$ of one or more output degrees of freedom.
- The perturbation amplitude sets the maximum recoverable range of the bifurcation diagram: the diagram can only be predicted up to amplitudes present in the data.
- Record for a sufficiently long time window that the system transits through the amplitude range of interest as it decays. This window is especially important close to the flutter boundary where recovery is slow.

---

## 4. Step-by-Step Algorithm

### Step 4.1: Mode Filtering (Large-Dimensional Systems)

In large systems, multiple modes contribute to the transient. Only the bifurcating mode exhibits CSD. Two approaches are available:

#### Option A: Bandpass Filtering (Simple)

1. Identify the frequency of the bifurcating mode from the **tail of the transient** at the parameter value closest to flutter. At late times all stable modes have decayed; the remaining oscillation frequency is $\omega_c$.
2. Apply a bandpass filter centred at $\omega_c$ to each transient signal. The filter width should be narrow enough to reject other modes but wide enough to pass the decaying envelope without distortion.
3. The filtered signal $q(t; \boldsymbol{\mu}_l)$ is used in subsequent steps.

This approach is unreliable when other modes have frequencies close to $\omega_c$.

#### Option B: Eigensystem Realisation Algorithm (ERA) Projection (Recommended)

1. Apply ERA to the transient responses at the parameter value $\boldsymbol{\mu}_l$ **closest to the flutter boundary**.
2. ERA constructs a state-space model from the Hankel matrix of measurements. Concretely:
   - Form the Hankel matrices $\mathbf{H}(k)$ from sampled output sequences $Y_k$:
     $$\mathbf{H}(k) = \begin{bmatrix} Y_k & Y_{k+1} & \cdots & Y_{k+\nu-1} \\ Y_{k+1} & Y_{k+2} & \cdots & Y_{k+\nu} \\ \vdots & & \ddots & \vdots \\ Y_{k+r-1} & Y_{k+r} & \cdots & Y_{k+r+\nu-2} \end{bmatrix}$$
   - Compute the SVD of $\mathbf{H}(0) = \mathbf{P}\mathbf{\Sigma}\mathbf{J}^T$. Retain the first $N$ singular values (model order).
   - Form truncated matrices $\bar{\mathbf{P}}, \bar{\mathbf{J}}, \bar{\mathbf{\Sigma}}$ from the first $N$ columns/values.
   - The discrete-time state-transition matrix and output matrix are:
     $$S = \bar{\mathbf{\Sigma}}^{-1/2}\bar{\mathbf{P}}^T\mathbf{H}(1)\bar{\mathbf{J}}\bar{\mathbf{\Sigma}}^{-1/2}, \quad C = \mathbf{E}_n^T\bar{\mathbf{P}}\bar{\mathbf{\Sigma}}^{1/2}$$
   - Diagonalise $S$ to obtain eigenvalues $\tilde{s}_i$ and eigenvectors. The continuous-time eigenvalue is $\eta_i = \frac{1}{\Delta t}\ln(\tilde{s}_i) = \sigma_i \pm j\omega_d$.
3. Identify the **bifurcating mode** as the eigenvalue with the most positive real part (least negative damping).
4. Extract the corresponding right eigenvector $\phi^{(c)}$ (complex, $N \times 1$) and left eigenvector $\psi^{(c)}$.
5. Form the $n_{\mathrm{out}} \times 2$ real matrices:
   $$\boldsymbol{\Phi} = \begin{bmatrix} \phi^{R} & \phi^{I} \end{bmatrix}, \quad \boldsymbol{\Psi} = \begin{bmatrix} \psi^{R} & \psi^{I} \end{bmatrix}$$

   **Biorthogonal left vector in physical output space.** To guarantee exact cancellation of the other ERA modes from the projected coordinate, $\psi^{(c)}$ must be biorthogonal to the physical output eigenvectors of all retained modes, not merely to their ERA state-space representations. Concretely: collect the first-block physical eigenvectors of every positive-imaginary ERA mode into $\boldsymbol{\Phi}_{\mathrm{all}} \in \mathbb{C}^{n_{\mathrm{out}} \times N_{\mathrm{modes}}}$, then form $\boldsymbol{\Psi}_{\mathrm{all}} = (\boldsymbol{\Phi}_{\mathrm{all}}^H)^+$ (Moore–Penrose pseudoinverse; equals $(\boldsymbol{\Phi}_{\mathrm{all}}^H)^{-1}$ when $N_{\mathrm{modes}} = n_{\mathrm{out}}$). Extract $\psi^{(c)}$ as the column of $\boldsymbol{\Psi}_{\mathrm{all}}$ corresponding to the bifurcating mode. Simple truncation of the ERA state-space left eigenvector does not preserve this biorthogonality.

6. The ERA is performed only once (at the nearest-to-flutter parameter value) and the resulting eigenvectors are reused for all other parameter samples, justified by the assumption that the centre space varies slowly with $\boldsymbol{\mu}$.
7. Project the transient response at each $\boldsymbol{\mu}_l$ onto the bifurcating mode using the complex inner-product form:
   $$q(t; \boldsymbol{\mu}_l) = \mathrm{Re}\!\left[\frac{(\psi^{(c)})^H\, \mathbf{y}(t;\boldsymbol{\mu}_l)}{(\psi^{(c)})^H\, \phi^{(c)}}\right]$$

   This scalar coordinate oscillates at $\omega_c$ and decays at rate $g_c$. The denominator $(\psi^{(c)})^H\phi^{(c)}$ has magnitude $\approx 1$ and is well-conditioned. The equivalent real $2\times 2$ form $(\boldsymbol{\Psi}^T\boldsymbol{\Phi})^{-1}\boldsymbol{\Psi}^T\mathbf{y}$ becomes catastrophically ill-conditioned when $\psi^{(c)}$ and $\phi^{(c)}$ are nearly orthogonal in real space — which occurs in lightly damped systems where the eigenvectors are predominantly imaginary — even though complex biorthogonality is maintained.

8. The projected signal $q(t; \boldsymbol{\mu}_l)$ replaces the raw output in subsequent steps.

---

### Step 4.2: Phase Fixing

Because the system oscillates during recovery, consecutive data points correspond to different oscillation phases and cannot be directly compared. Phase fixing selects a consistent phase.

**Standard approach:** Extract the **local maxima** of $|q(t; \boldsymbol{\mu}_l)|$. Denote the peak times and values as $(t_k, r_k)$. These correspond to the phase where $\dot{q} = 0$, $q > 0$. The local minima may be used separately to forecast the lower branch of asymmetric bifurcation diagrams.

**Note:** Discard the initial portion of the transient (e.g.\ the first oscillation period or one second) to ensure the data lies on the inertial manifold after initial transients from non-bifurcating modes have decayed.

---

### Step 4.3: Estimating the Recovery Rate

Two methods are available depending on the number of available peaks.

#### Option A: Finite Differences (High-Frequency Oscillations, Many Peaks)

For each consecutive pair of peaks $(t_k, r_k)$ and $(t_{k+1}, r_{k+1})$:

$$\lambda_k = \frac{\ln r_{k+1} - \ln r_k}{t_{k+1} - t_k}$$

This gives a cloud of $(\lambda_k, r_k)$ pairs for each $\boldsymbol{\mu}_l$.

A centred finite-difference scheme using the amplitudes before and after can also be used:

$$\lambda(t, \boldsymbol{\mu}_l) \approx \frac{\ln r_+ - \ln r_-}{2\Delta t}$$

where $r_\pm$ are amplitudes at $t \pm \Delta t$, applied only at peak times.

#### Option B: Nonlinear Optimisation (Low-Frequency Oscillations, Few Peaks)

When only a small number of peaks are available per transient, model the recovery rate at fixed $\boldsymbol{\mu}_l$ as a polynomial in $r$:

$$\lambda(r; \boldsymbol{\mu}_l) = \lambda_0 + \lambda_1 r + \lambda_2 r^2 + \cdots + \lambda_p r^p$$

This implies the amplitude ODE:

$$\dot{r} = r(\lambda_0 + \lambda_1 r + \cdots + \lambda_p r^p)$$

**Procedure:**
1. Choose polynomial order $p$ (default $p = 4$; determined by a convergence test — increase $p$ until results stabilise).
2. Integrate the ODE forward in time from the initial perturbation amplitude $r_0$ for a candidate set of coefficients $\{\lambda_i\}$.
3. Minimise the sum-of-squares mismatch between the integrated envelope $r(t_k)$ and the observed peak values $r_k^{\mathrm{obs}}$ using nonlinear least squares (e.g.\ Levenberg–Marquardt):
   $$\min_{\lambda_0,\ldots,\lambda_p} \sum_k \left(r(t_k; \lambda_0,\ldots,\lambda_p) - r_k^{\mathrm{obs}}\right)^2$$
4. The fitted polynomial $\lambda(r; \boldsymbol{\mu}_l)$ gives recovery rate values at any amplitude for this parameter sample.

This approach recovers a smooth, continuous $\lambda(r)$ curve per parameter sample from as few as 3–4 peaks.

---

### Step 4.4: Building the Recovery Rate Dataset

After Steps 4.1–4.3, for each parameter sample $\boldsymbol{\mu}_l$ one has either:
- A discrete set of $(\lambda_k, r_k)$ pairs (Option A), or
- A smooth polynomial $\lambda(r; \boldsymbol{\mu}_l)$ evaluated at any desired $r$ (Option B).

In either case, select a grid of amplitude values $\{\tilde{r}_1, \tilde{r}_2, \ldots, \tilde{r}_{N_r}\}$ spanning the range of interest. For each $\tilde{r}_k$ and each $\boldsymbol{\mu}_l$, interpolate or evaluate to obtain:

$$\lambda_{kl} = \lambda(\tilde{r}_k, \boldsymbol{\mu}_l)$$

This gives an $N_r \times N_s$ matrix of recovery rates indexed by amplitude and parameter sample.

**Important:** Only amplitude values below the smallest initial perturbation across all simulations can be predicted. Restrict $\tilde{r}_k$ accordingly.

---

### Step 4.5: Fitting the Recovery Rate Surface

At each fixed amplitude $\tilde{r}_k$, fit a polynomial in $\boldsymbol{\mu}$ to the data $\{(\boldsymbol{\mu}_l, \lambda_{kl})\}_{l=1}^{N_s}$.

The fitting polynomial (working in raw $\boldsymbol{\mu}$ coordinates, absorbing the shift from $\tilde{\boldsymbol{\mu}}_c$) is:

$$\lambda(\tilde{r}_k, \boldsymbol{\mu}) \approx a_0(\tilde{r}_k) + \mathbf{a}_1(\tilde{r}_k)^T \boldsymbol{\mu} + \boldsymbol{\mu}^T \mathbf{a}_2(\tilde{r}_k)\, \boldsymbol{\mu}$$

where:
- $a_0(\tilde{r}_k) \in \mathbb{R}$: zeroth-order coefficient.
- $\mathbf{a}_1(\tilde{r}_k) \in \mathbb{R}^{N_p}$: first-order coefficients (one per parameter).
- $\mathbf{a}_2(\tilde{r}_k) \in \mathbb{R}^{N_p \times N_p}$: second-order coefficient matrix (symmetric).

**Explicit polynomial terms for $N_p = 2$, $\boldsymbol{\mu} = (U, \sigma)$:**

*First-order (bi-linear, $N_s \geq 3$):*
$$\lambda \approx a_0 + a_1^{(1)} U + a_1^{(2)} \sigma$$

*Second-order (bi-quadratic, $N_s \geq 6$):*
$$\lambda \approx a_0 + a_1^{(1)} U + a_1^{(2)} \sigma + a_2^{(11)} U^2 + a_2^{(12)} U\sigma + a_2^{(22)} \sigma^2$$

**Fitting procedure:**

For each amplitude $\tilde{r}_k$:
1. Write one scalar equation of the form above for each $\boldsymbol{\mu}_l$.
2. Stack into a linear system $\mathbf{A}_k \mathbf{c}_k = \boldsymbol{\lambda}_k$ where:
   - $\mathbf{A}_k$ is $N_s \times N_c$ (design matrix of monomial values at each sample point), $N_c$ = number of polynomial coefficients.
   - $\mathbf{c}_k$ is the $N_c \times 1$ vector of unknown coefficients.
   - $\boldsymbol{\lambda}_k = [\lambda_{k1}, \ldots, \lambda_{kN_s}]^T$.
3. Solve by least squares: $\mathbf{c}_k = (\mathbf{A}_k^T \mathbf{A}_k)^{-1} \mathbf{A}_k^T \boldsymbol{\lambda}_k$ (or via QR factorisation for numerical stability).
4. Store the coefficient vector $\mathbf{c}_k$ for each amplitude $\tilde{r}_k$.

**Practical notes:**
- Normalise parameters before fitting: subtract the mean and divide by the range to condition the design matrix.
- Start with first-order expansion; use a convergence check (compare predicted flutter speed under first- vs.\ second-order) to decide whether higher order is needed.
- Second-order fitting is more sensitive to noise and numerical errors; use additional transient samples to stabilise it.

---

### Step 4.6: Extracting the Bifurcation Diagram

For each amplitude $\tilde{r}_k$, find the control parameter values where $\lambda(\tilde{r}_k, \boldsymbol{\mu}) = 0$.

#### Single parameter ($N_p = 1$)

With $\boldsymbol{\mu} = \mu$ (scalar), the fitted polynomial is:

$$\lambda(\tilde{r}_k, \mu) = a_0 + a_1 \mu + a_2 \mu^2 = 0$$

Solve the quadratic (or linear if second-order terms are dropped) analytically:

$$\tilde{\mu}_k = \frac{-a_1 \pm \sqrt{a_1^2 - 4a_2 a_0}}{2a_2}$$

Choose the root closest to the pre-flutter samples and physically meaningful (i.e.\ $\tilde{\mu}_k > \mu_l$ for flutter onset beyond the measured speeds).

#### Multiple parameters ($N_p \geq 2$)

Parameterise the search direction in parameter space. Write the bifurcation diagram point as:

$$\tilde{\boldsymbol{\mu}} = \boldsymbol{\mu}^\star + \tilde{\delta}\, \bar{\boldsymbol{\mu}}$$

where $\boldsymbol{\mu}^\star$ is a chosen reference point, $\bar{\boldsymbol{\mu}}$ is a unit direction vector, and $\tilde{\delta}$ is the scalar unknown.

Substituting into the zero-recovery-rate condition gives a **scalar quadratic in $\tilde{\delta}$**:

$$a_0 + \mathbf{a}_1^T(\boldsymbol{\mu}^\star + \tilde{\delta}\bar{\boldsymbol{\mu}}) + (\boldsymbol{\mu}^\star + \tilde{\delta}\bar{\boldsymbol{\mu}})^T \mathbf{a}_2 (\boldsymbol{\mu}^\star + \tilde{\delta}\bar{\boldsymbol{\mu}}) = 0$$

Expanding:

$$\left(\bar{\boldsymbol{\mu}}^T \mathbf{a}_2 \bar{\boldsymbol{\mu}}\right)\tilde{\delta}^2 + \left(\mathbf{a}_1^T\bar{\boldsymbol{\mu}} + 2(\boldsymbol{\mu}^\star)^T\mathbf{a}_2\bar{\boldsymbol{\mu}}\right)\tilde{\delta} + \left(a_0 + \mathbf{a}_1^T\boldsymbol{\mu}^\star + (\boldsymbol{\mu}^\star)^T\mathbf{a}_2\boldsymbol{\mu}^\star\right) = 0$$

Solve analytically for $\tilde{\delta}$; the bifurcation diagram point is $\tilde{\boldsymbol{\mu}} = \boldsymbol{\mu}^\star + \tilde{\delta}\bar{\boldsymbol{\mu}}$.

**Tracing the flutter boundary contour:**
- Fix $\bar{\boldsymbol{\mu}}$ and sweep $\boldsymbol{\mu}^\star$ to obtain parallel sections of the boundary.
- Fix $\boldsymbol{\mu}^\star$ and sweep $\bar{\boldsymbol{\mu}}$ (varying the search angle) to find the most critical direction.

Repeat for all $\tilde{r}_k$ to build the full multi-dimensional bifurcation diagram in $(\boldsymbol{\mu}, r)$ space.

---

### Step 4.7: Identifying Bifurcation Type and LCO Stability

- **Flutter boundary:** Set $r = 0$ (equivalently, $a_0 = 0$ from the zero-$r$ limit of the fitting polynomial). Solve for $\boldsymbol{\mu}_c$.
- **Bifurcation type:** Two approaches are available.

  *Primary (robust) — majority-vote position comparison:* For each grid amplitude $\tilde{r}_k$, compute the bifurcation diagram point $\tilde{\boldsymbol{\mu}}_k$ and its signed displacement from the flutter boundary along the search direction $\hat{\boldsymbol{\mu}}$:
  $$\delta_k = (\tilde{\boldsymbol{\mu}}_k - \boldsymbol{\mu}_c)\cdot\hat{\boldsymbol{\mu}}$$
  Classify by majority vote over all $k$: if most $\delta_k > 0$ the branch lies above the flutter boundary (supercritical); if most $\delta_k < 0$ the branch lies below (subcritical). Using all diagram points rather than just small-$r$ ones is essential: at small $r$ the LCO shift $\sim r^2$ is comparable to the surface-fitting residual and the per-point comparison is unreliable; at large $r$ the branch is unambiguously separated from the flutter boundary.

  *Fallback — $\partial\lambda/\partial r$ sign:* If the flutter boundary is unavailable, classify from the majority vote on the sign of $\partial\lambda/\partial r|_{\lambda=0}$: negative slope (supercritical); positive slope (subcritical). This is less robust than position comparison for the same reason.

- **LCO stability:** The slope $\partial\lambda/\partial r|_{\lambda=0}$ at a bifurcation diagram point determines stability of that LCO solution. Negative slope = stable LCO.
- **LCO frequency:** The oscillation frequency of the post-flutter LCO is not predicted by this framework. The linear flutter frequency $\omega_c$ (imaginary part of the bifurcating mode eigenvalue at the flutter boundary) is available as a by-product of the ERA or eigenvalue analysis, but frequency-amplitude variation in the post-flutter regime requires separate analysis (e.g.\ harmonic balance, multiple scales).

---

## 5. Complete Algorithm Summary (Pseudocode)

```
INPUT:  Parameter samples {μ_l}, l = 1..Ns
        Free-decay responses y(t; μ_l) at each sample

--- PREPROCESSING ---
ERA on response at μ closest to flutter
  → bifurcating mode eigenvectors Φ, Ψ
  → flutter frequency ω_c

For each l = 1..Ns:
  Project y(t; μ_l) onto bifurcating mode → q(t; μ_l)
  Discard initial transient (≈ 1 oscillation period)
  Extract local maxima → peaks (t_k, r_k)

--- RECOVERY RATE ESTIMATION ---
For each l = 1..Ns:
  If many peaks (>~10):
    λ_kl = (ln r_{k+1} - ln r_k) / (t_{k+1} - t_k)   for each k
  Else (few peaks):
    Fit ODE: r_dot = r(λ_0 + λ_1*r + ... + λ_p*r^p) to peaks by nonlinear least squares
    Evaluate λ(r; μ_l) at desired amplitudes

--- SURFACE FITTING ---
Select amplitude grid {r̃_1, ..., r̃_Nr}
For each amplitude r̃_k:
  Assemble design matrix A_k  (rows = samples μ_l, cols = polynomial monomials)
  Assemble RHS vector λ_k = [λ(r̃_k, μ_1), ..., λ(r̃_k, μ_Ns)]^T
  Solve least-squares: c_k = argmin ||A_k c - λ_k||²
  Store coefficients c_k = {a_0(r̃_k), a_1(r̃_k), a_2(r̃_k)}

--- BIFURCATION DIAGRAM ---
Flutter boundary:
  At r̃ = 0: solve a_0(0) + a_1(0)^T μ + μ^T a_2(0) μ = 0 for μ_c
  (Trace contour in parameter space by sweeping search directions)

For each r̃_k > 0:
  For each search direction μ_bar and reference point μ_star:
    Solve scalar quadratic in δ̃:
      (μ_bar^T a_2 μ_bar) δ̃² + (a_1^T μ_bar + 2 μ_star^T a_2 μ_bar) δ̃
        + (a_0 + a_1^T μ_star + μ_star^T a_2 μ_star) = 0
    Bifurcation point: μ̃ = μ_star + δ̃ * μ_bar

OUTPUT: Flutter boundary in μ-space
        Bifurcation diagram: (μ̃, r̃) pairs
        Bifurcation type (super/subcritical): majority vote over all diagram points —
          LCO branch above flutter boundary (δ_k > 0) → supercritical;
          below (δ_k < 0) → subcritical; fallback to sign of ∂λ/∂r if boundary unavailable
```

---

## 6. Extension: State Velocity Method (Riso, Cesnik & Epureanu, AIAA J 2022)

This variant **replaces transient simulations with direct evaluation of the state velocity**, making it suitable for model-based design applications. It requires a time-domain state-space model $\dot{\mathbf{y}} = \mathbf{f}(p, \mathbf{y})$.

### 6.1 Formulation

For each control parameter value $p_i$ and amplitude $r_k$:

1. Find the equilibrium state $\mathbf{y}_{e_i}$ satisfying $\mathbf{f}(p_i, \mathbf{y}_{e_i}) = \mathbf{0}$.
2. Compute the Jacobian $\mathbf{A}_i = \partial\mathbf{f}/\partial\mathbf{y}|_{p_i, \mathbf{y}_{e_i}}$ and find the eigenvectors $\boldsymbol{\phi}_{ic}$, $\boldsymbol{\psi}_{ic}$ of the bifurcating mode (most positive real part of eigenvalue).
3. Form the $N \times 2$ matrices $\boldsymbol{\Phi}_{ic} = [\boldsymbol{\phi}^R_{ic},\, \boldsymbol{\phi}^I_{ic}]$ and $\boldsymbol{\Psi}_{ic} = [\boldsymbol{\psi}^R_{ic},\, \boldsymbol{\psi}^I_{ic}]$.
4. Reduce the $N \times N$ state matrix to a $2 \times 2$ critical-mode matrix:
   $$\tilde{\mathbf{A}}_{ic} = (\boldsymbol{\Psi}_{ic}^T \boldsymbol{\Phi}_{ic})^{-1} \boldsymbol{\Psi}_{ic}^T \mathbf{A}_i \boldsymbol{\Phi}_{ic} = \begin{bmatrix} g_{ic} & \omega_{ic} \\ -\omega_{ic} & g_{ic} \end{bmatrix}$$
5. Construct the perturbed state at amplitude $r$ along the bifurcating mode at phase $\theta = 0$:
   $$\mathbf{y}_{0_{ic}}(r) = \mathbf{y}_{e_i} + \boldsymbol{\Phi}_{ic}\tilde{\boldsymbol{\phi}}^R_{ic}\, r$$
   where $\tilde{\boldsymbol{\phi}}^R_{ic} = \{1, 0\}^T$ (real part of reduced eigenvector).
6. Evaluate the full state velocity at this perturbed state:
   $$\boldsymbol{\Phi}_{ic}\tilde{\boldsymbol{\phi}}^R_{ic}\dot{r} = \mathbf{f}(p_i,\, \mathbf{y}_{0_{ic}}(r))$$
7. Project onto the reduced left eigenvectors to obtain the scalar amplitude rate:
   $$\dot{r} = (\tilde{\boldsymbol{\psi}}^R_{ic})^T(\boldsymbol{\Psi}_{ic}^T\boldsymbol{\Phi}_{ic})^{-1}\boldsymbol{\Psi}_{ic}^T\, \mathbf{f}(\boldsymbol{\Phi}_{ic}\tilde{\boldsymbol{\phi}}^R_{ic}\, r) \equiv \tilde{f}_{ic}(r)$$
8. The recovery rate is then:
   $$\lambda_{ic}(r) \approx \frac{\tilde{f}_{ic}(r)}{r} = g_{ic} + (\tilde{\boldsymbol{\psi}}^R_{ic})^T(\boldsymbol{\Psi}_{ic}^T\boldsymbol{\Phi}_{ic})^{-1}\boldsymbol{\Psi}_{ic}^T \frac{\mathbf{f}_{nl}(\boldsymbol{\Phi}_{ic}\tilde{\boldsymbol{\phi}}^R_{ic}\, r)}{r}$$
   where $\mathbf{f}_{nl}$ is the nonlinear part of $\mathbf{f}$ (i.e.\ $\mathbf{f} - \mathbf{A}_i\mathbf{y}$).

Steps 1–8 are repeated at each $(p_i, r_k)$ without any time integration.

### 6.2 Fitting and Diagram Extraction

Identical to Steps 4.5–4.6 of the primary method. Recovery rates $\lambda_{ik}$ are fitted as a polynomial in $p$ at each amplitude, then zero-crossing gives the bifurcation diagram.

### 6.3 Accuracy and Limitations

- Provides exact recovery rates for polynomial systems (e.g.\ cubic/quintic stiffness); approximate for more complex nonlinearities.
- Captures bifurcation type (super/subcritical) and approximate bifurcation diagram amplitude trend; quantitative accuracy is lower than the transient-data method.
- Requires only eigenvalue analyses and function evaluations — no nonlinear time integration.
- Inaccurate if nonlinear normal modes vary strongly with amplitude (eigenvectors at equilibrium become a poor basis).
- Not applicable when a model is unavailable.

---

## 7. Extension: MMS-Hybrid Single-Trajectory Method (García Pérez et al., 2023)

A further extension by García Pérez, Ghadami, Sanches, Michon & Epureanu combines bifurcation forecasting with the **method of multiple scales (MMS)**. Rather than requiring multiple fixed-parameter transients, this approach uses a single trajectory in which the bifurcation parameter varies continuously during the recovery.

The MMS provides an asymptotic normal form capturing CSD as $\boldsymbol{\mu}$ changes during the recovery. This removes the requirement for separate fixed-parameter runs. Accuracy is improved at large amplitudes compared to the standard method, but the approach requires knowledge of the system's nonlinear structure to obtain the normal form.

---

## 8. Practical Implementation Notes

### Parameter Grid Design

- Use LHS to generate parameter samples in a box centred on the expected flutter region.
- For first-order fitting with $N_p = 2$: 3 samples minimum; 5–6 recommended for noise robustness.
- For second-order fitting with $N_p = 2$: 6 samples minimum; 8–10 recommended.
- Samples should span variations in each parameter independently and in combination.

### Amplitude Grid Selection

- Select $N_r = 20$–50 amplitude values spanning from near-zero to the maximum perturbation amplitude.
- Use finer spacing near $r = 0$ to resolve the flutter boundary accurately.
- Restrict to amplitudes below the initial perturbation of every simulation used.
- **Common-overlap region.** When curves from different parameter samples cover different amplitude ranges (e.g.\ samples close to the flutter boundary decay slowly and have a high noise floor at small $r$), restrict the grid to the intersection of all curve ranges:
  $$r_{\min} = \max_l\!\left(\min_k r_{lk}\right)\!\times 1.05, \qquad r_{\max} = \min_l\!\left(\max_k r_{lk}\right)\!\times 0.95$$
  Using amplitudes below the lowest curve's minimum forces extrapolation that produces unreliable zero-crossing estimates and spurious bifurcation-type classifications, especially at small $r$ where the surface fit residual dominates the nonlinear correction $\sim r^2$.

### Polynomial Fitting Stability

- Normalise each parameter: $\mu_i \leftarrow (\mu_i - \bar{\mu}_i)/\Delta\mu_i$ before assembling the design matrix.
- Check condition number of $\mathbf{A}_k^T\mathbf{A}_k$; if ill-conditioned, add Tikhonov regularisation: $\mathbf{c}_k = (\mathbf{A}_k^T\mathbf{A}_k + \varepsilon\mathbf{I})^{-1}\mathbf{A}_k^T\boldsymbol{\lambda}_k$.
- The fitting coefficients $a_0(\tilde{r})$, $a_1^{(i)}(\tilde{r})$, $a_2^{(ij)}(\tilde{r})$ will be smooth functions of amplitude; if they are noisy, apply mild smoothing before the extrapolation step.

### Root-Finding

- The scalar quadratic in $\tilde{\delta}$ has a closed-form solution.
- If the quadratic has no real roots at some amplitudes, the bifurcation diagram cannot be predicted there (either out of extrapolation range or the quadratic approximation has broken down).
- For the flutter boundary ($r = 0$), use only the linear term of the amplitude polynomial to avoid noise amplification.

### Noise and Robustness

- With 2–5% measurement noise, increase $N_s$ (more transient samples) rather than increasing polynomial order.
- The nonlinear optimisation approach (Step 4.3 Option B) is inherently more robust to noise than finite differences.
- ERA-based modal filtering is particularly effective at suppressing measurement noise as a side effect.

### Convergence Check

1. Run with first-order Taylor expansion; record predicted flutter speed $\tilde{\mu}_c^{(1)}$.
2. Run with second-order; record $\tilde{\mu}_c^{(2)}$.
3. If $|\tilde{\mu}_c^{(2)} - \tilde{\mu}_c^{(1)}|/|\tilde{\mu}_c^{(1)}| < \varepsilon_{\mathrm{tol}}$ (e.g.\ 1%), first order is sufficient.
4. Similarly compare bifurcation diagram shapes; if they agree at small amplitudes, first order is adequate.

---

## 9. Key Equations Reference

| Quantity | Expression |
|---|---|
| Recovery rate definition | $\lambda = \dot{r}/r = \mathrm{d}(\ln r)/\mathrm{d}t$ |
| Finite-difference estimate | $\lambda_k = (\ln r_{k+1} - \ln r_k)/(t_{k+1} - t_k)$ |
| Nonlinear ODE for envelope | $\dot{r} = r\sum_{n=0}^p \lambda_n r^n$ |
| Taylor expansion (general) | $\lambda \approx \alpha_0 + \boldsymbol{\alpha}_1^T(\boldsymbol{\mu}-\tilde{\boldsymbol{\mu}}_c) + (\boldsymbol{\mu}-\tilde{\boldsymbol{\mu}}_c)^T\boldsymbol{\alpha}_2(\boldsymbol{\mu}-\tilde{\boldsymbol{\mu}}_c)$ |
| Fitting polynomial (raw $\boldsymbol{\mu}$) | $\lambda \approx a_0 + \mathbf{a}_1^T\boldsymbol{\mu} + \boldsymbol{\mu}^T\mathbf{a}_2\boldsymbol{\mu}$ |
| Bifurcation diagram condition | $\lambda(\tilde{\boldsymbol{\mu}}, \tilde{r}) = 0$ |
| Parametric search | $\tilde{\boldsymbol{\mu}} = \boldsymbol{\mu}^\star + \tilde{\delta}\bar{\boldsymbol{\mu}}$ |
| Quadratic in $\tilde{\delta}$ | $(\bar{\boldsymbol{\mu}}^T\mathbf{a}_2\bar{\boldsymbol{\mu}})\tilde{\delta}^2 + (\mathbf{a}_1^T\bar{\boldsymbol{\mu}} + 2(\boldsymbol{\mu}^\star)^T\mathbf{a}_2\bar{\boldsymbol{\mu}})\tilde{\delta} + (a_0 + \mathbf{a}_1^T\boldsymbol{\mu}^\star + (\boldsymbol{\mu}^\star)^T\mathbf{a}_2\boldsymbol{\mu}^\star) = 0$ |
| Minimum samples (1st order) | $N_s \geq 1 + N_p$ |
| Minimum samples (2nd order) | $N_s \geq 1 + N_p + N_p(N_p+1)/2$ |
| State velocity recovery rate | $\lambda_{ic}(r) = g_{ic} + (\tilde{\boldsymbol{\psi}}^R)^T(\boldsymbol{\Psi}^T\boldsymbol{\Phi})^{-1}\boldsymbol{\Psi}^T\mathbf{f}_{nl}(\boldsymbol{\Phi}\tilde{\boldsymbol{\phi}}^R r)/r$ |

---

## 10. Primary References

| Reference | DOI | Contribution |
|---|---|---|
| Lim & Epureanu (2011) *Phys. Rev. E* | `10.1103/PhysRevE.83.016203` | Original bifurcation forecasting for non-oscillatory systems |
| Ghadami & Epureanu (2016) *J. Comput. Nonlinear Dyn.* | `10.1115/1.4033920` | Extension to large-dimensional oscillatory systems; ERA modal filtering; nonlinear optimisation for recovery rates |
| Ghadami & Epureanu (2017) *Nonlinear Dyn.* | `10.1007/s11071-016-3250-y` | Centre-space projection; slow-oscillatory systems |
| Ghadami, Cesnik & Epureanu (2018) *J. Fluids Struct.* | `10.1016/j.jfluidstructs.2017.09.005` | Application to geometrically nonlinear high-aspect-ratio wing |
| Ghadami & Epureanu (2018) *Int. J. Non-Linear Mech.* | `10.1016/j.ijnonlinmec.2018.02.008` | Forecasting critical points and post-critical LCOs |
| Riso, Ghadami, Cesnik & Epureanu (2020) *AIAA J.* | `10.2514/1.J059024` | Single-parameter application to geometrically nonlinear wings |
| Riso, Cesnik & Epureanu (2021) *J. Fluids Struct.* | `10.1016/j.jfluidstructs.2020.103201` | **Multi-parameter extension**; LHS sampling; bi-linear/quadratic fitting |
| Riso, Cesnik & Epureanu (2022) *AIAA J.* | `10.2514/1.J061860` | State velocity method; model-based recovery rate without transients |
| Gali, Goehmann & Riso (2023) *J. Fluids Struct.* | `10.1016/j.jfluidstructs.2023.103948` | Application to whirl flutter |
| García Pérez et al. (2023) *Nonlinear Dyn.* | `10.1007/s11071-023-08502-x` | MMS-hybrid single-trajectory method |
