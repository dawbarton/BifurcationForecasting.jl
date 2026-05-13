# Test Case 1: Analytical 1-DOF Bifurcation Forecasting

*Based on Riso, Cesnik & Epureanu, Journal of Fluids and Structures 101 (2021) 103201, Section 3.*

---

## 1. Purpose

This test case verifies the core bifurcation forecasting pipeline — transient simulation, recovery rate
estimation, polynomial surface fitting, and zero-crossing extrapolation — against exact analytical
solutions. Because the recovery rate is a known polynomial in the control parameters, a correct
implementation must reproduce the bifurcation diagram to within numerical integration error.

---

## 2. Governing Equations

Two scalar amplitude ODEs are studied. Both have equilibrium at $r = 0$ and bifurcate at $\mu = 0$.

### Supercritical system

$$\dot{r} = r(\alpha\mu - \beta r^2) \tag{1}$$

A stable LCO branch is born at $\mu = 0$ and grows for $\mu > 0$.

### Subcritical system

$$\dot{r} = r(\alpha\mu + \beta r^2 - \gamma r^4) \tag{2}$$

An unstable LCO branch exists for $\mu < 0$; a fold point occurs at $\mu = -\beta^2/(4\alpha\gamma)$.

### Fixed parameters

$$\alpha = 1, \quad \gamma = 1$$

The two **varying** (control) parameters are $\mu$ and $\beta$. Both are treated as free in the
two-parameter forecasting demonstration.

---

## 3. Analytical Solutions

### Recovery rates

From $\lambda = \dot{r}/r$:

$$\lambda_{\mathrm{sup}}(\mu, \beta, r) = \alpha\mu - \beta r^2 \tag{3}$$

$$\lambda_{\mathrm{sub}}(\mu, \beta, r) = \alpha\mu + \beta r^2 - \gamma r^4 \tag{4}$$

### Bifurcation diagram (zero of $\lambda$)

**Supercritical stable branch** ($\mu \geq 0$, $\beta > 0$, $\gamma = 0$):

$$\tilde{r}_{\mathrm{sup}}(\mu, \beta) = \sqrt{\frac{\alpha\mu}{\beta}} \tag{5}$$

**Subcritical stable branch** ($\mu \geq -\beta^2/(4\alpha\gamma)$, $\beta > 0$, $\gamma > 0$):

$$\tilde{r}_{\mathrm{sub,1}}(\mu, \beta, \gamma) = \sqrt{\frac{\beta + \sqrt{\beta^2 + 4\alpha\gamma\mu}}{2\gamma}} \tag{6}$$

**Subcritical unstable branch** ($\mu \in [-\beta^2/(4\alpha\gamma),\, 0]$):

$$\tilde{r}_{\mathrm{sub,2}}(\mu, \beta, \gamma) = \sqrt{\frac{\beta - \sqrt{\beta^2 + 4\alpha\gamma\mu}}{2\gamma}} \tag{7}$$

**Flutter boundary** (both systems): $\mu_c = 0$ for all $\beta$.

---

## 4. Simulated Transients

### Initial condition

$$r_0 = 1.5 \quad \text{for all simulations}$$

### ODE integration

Integrate the chosen ODE (Eq. 1 or Eq. 2) forward in time from $r(0) = r_0$ using a standard
fixed-step or adaptive solver (e.g. RK4 with step $\Delta t \leq 0.01$). Because these are scalar
non-oscillatory ODEs, no phase fixing or modal filtering is required — the signal $r(t)$ is the
amplitude envelope directly.

Record $r(t)$ at uniform time steps until $r < 10^{-4}$ or until $t = t_{\max}$ (choose
$t_{\max} = 50$ for samples far from flutter; up to $t_{\max} = 200$ for samples close to $\mu = 0$).

### Parameter samples

The following nine samples reproduce Table 1 of Riso JFS 2021:

| Simulation | $\mu$ | $\beta$ | Purpose |
|---|---|---|---|
| 1 | −0.50 | 1.00 | Single-parameter fitting ($\beta$ fixed) |
| 2 | −0.40 | 1.00 | Single-parameter fitting |
| 3 | −0.30 | 1.00 | Single-parameter fitting |
| 4 | −0.48 | 1.05 | Two-parameter fitting |
| 5 | −0.36 | 1.13 | Two-parameter fitting |
| 6 | −0.32 | 0.97 | Two-parameter fitting |
| 7 | −0.39 | 0.84 | Two-parameter fitting |
| 8 | −0.43 | 0.92 | Two-parameter fitting |
| 9 | −0.43 | 1.13 | Two-parameter fitting |

All samples lie in the pre-flutter regime ($\mu < 0$).

---

## 5. Recovery Rate Estimation

Because these systems are non-oscillatory, no phase fixing is needed. The recovery rate is estimated
directly from the time series using finite differences.

For consecutive time steps $(t_k, r_k)$ and $(t_{k+1}, r_{k+1})$:

$$\lambda_k = \frac{\ln r_{k+1} - \ln r_k}{t_{k+1} - t_k} \tag{8}$$

This gives a set of $(\lambda_k, r_k)$ pairs for each simulation $l$.

**Verification step A:** Compare the estimated $\lambda_k$ values against the analytical expressions
(Eqs. 3–4) evaluated at $(r_k, \mu_l, \beta_l)$. The maximum absolute error should be of order
$O(\Delta t)$ (first-order finite differences) or $O(\Delta t^2)$ with centred differences. For
$\Delta t = 0.01$ and a well-resolved trajectory, errors below $10^{-3}$ are expected.

---

## 6. Amplitude Grid

Select $N_r = 50$ amplitude values uniformly spaced in $[r_{\min}, r_{\max}]$:

$$r_{\min} = 0.01, \quad r_{\max} = 1.45$$

The upper limit is set slightly below $r_0 = 1.5$ to avoid edge effects at the start of the
transient. The lower limit avoids division by zero in recovery rate estimates.

For each $\tilde{r}_k$ and each simulation $l$, interpolate the $(r, \lambda)$ data from Step 5
to obtain $\lambda_{kl} = \lambda(\tilde{r}_k, \mu_l, \beta_l)$. Use linear interpolation in $r$.

---

## 7. Polynomial Surface Fitting

### 7.1 Single-parameter case (simulations 1–3, $\beta = 1$ fixed)

At each amplitude $\tilde{r}_k$, fit a quadratic polynomial in $\mu$ using the three data points
$\{(\mu_l, \lambda_{kl})\}_{l=1}^{3}$:

$$\lambda(\tilde{r}_k, \mu) \approx a_0(\tilde{r}_k) + a_1(\tilde{r}_k)\,\mu + a_2(\tilde{r}_k)\,\mu^2 \tag{9}$$

This requires solving the $3 \times 3$ Vandermonde system:

$$\begin{bmatrix} 1 & \mu_1 & \mu_1^2 \\ 1 & \mu_2 & \mu_2^2 \\ 1 & \mu_3 & \mu_3^2 \end{bmatrix} \begin{bmatrix} a_0 \\ a_1 \\ a_2 \end{bmatrix} = \begin{bmatrix} \lambda_{k1} \\ \lambda_{k2} \\ \lambda_{k3} \end{bmatrix} \tag{10}$$

For the analytical systems, the true recovery rate is linear in $\mu$ (from Eqs. 3–4), so the fitted
$a_2$ should be negligibly small (of order numerical noise).

A **linear** fit (dropping $a_2$) uses only two simulations (1 and 2, or 1 and 3) and gives:

$$a_0 = \frac{\mu_2 \lambda_{k1} - \mu_1 \lambda_{k2}}{\mu_2 - \mu_1}, \quad a_1 = \frac{\lambda_{k2} - \lambda_{k1}}{\mu_2 - \mu_1} \tag{11}$$

### 7.2 Two-parameter case (simulations 1–9)

The control-parameter vector is $\boldsymbol{\mu} = (\mu, \beta)^T$.

**First-order (bi-linear) fit** using simulations 4–6 (minimum 3 samples for $N_p = 2$, or
include samples 1–3 for 6 total):

$$\lambda(\tilde{r}_k, \mu, \beta) \approx a_0 + a_1^{(\mu)}\mu + a_1^{(\beta)}\beta \tag{12}$$

Design matrix ($N_s \times 3$) with one row per simulation:

$$\mathbf{A} = \begin{bmatrix} 1 & \mu_1 & \beta_1 \\ \vdots & \vdots & \vdots \\ 1 & \mu_{N_s} & \beta_{N_s} \end{bmatrix} \tag{13}$$

Solve by least squares: $\mathbf{c} = (\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T\boldsymbol{\lambda}_k$.

**Second-order (bi-quadratic) fit** using simulations 4–9 (minimum 6 samples):

$$\lambda \approx a_0 + a_1^{(\mu)}\mu + a_1^{(\beta)}\beta + a_2^{(\mu\mu)}\mu^2 + a_2^{(\mu\beta)}\mu\beta + a_2^{(\beta\beta)}\beta^2 \tag{14}$$

Design matrix ($N_s \times 6$):

$$\mathbf{A} = \begin{bmatrix} 1 & \mu_1 & \beta_1 & \mu_1^2 & \mu_1\beta_1 & \beta_1^2 \\ \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \end{bmatrix} \tag{15}$$

**Note on normalisation.** Before assembling $\mathbf{A}$, normalise each parameter:

$$\hat{\mu} = \frac{\mu - \bar{\mu}}{\sigma_\mu}, \quad \hat{\beta} = \frac{\beta - \bar{\beta}}{\sigma_\beta}$$

where $\bar{\mu}, \bar{\beta}$ are sample means and $\sigma_\mu, \sigma_\beta$ are sample standard
deviations. Solve for coefficients in normalised coordinates; transform back for evaluation.

**Verification step B:** For the analytical systems, the true recovery rate is exactly bi-linear in
$(\mu, \beta)$ (Eqs. 3–4). Therefore, the first-order fit should reproduce the recovery rate
**exactly** (up to floating-point and interpolation error). The second-order fit should give
identical results to the first-order fit, with the second-order coefficients being of order
machine epsilon.

---

## 8. Bifurcation Diagram Extraction

### 8.1 Single-parameter case

At each amplitude $\tilde{r}_k$, find the value $\tilde{\mu}_k$ where $\lambda(\tilde{r}_k, \tilde{\mu}_k) = 0$.

From the quadratic fit (Eq. 9):

$$\tilde{\mu}_k = \frac{-a_1 - \sqrt{a_1^2 - 4a_2 a_0}}{2a_2} \tag{16}$$

(choose the root that lies above all sample $\mu$ values, i.e.\ the root closer to zero from below).

If using the linear fit ($a_2 = 0$):

$$\tilde{\mu}_k = -\frac{a_0}{a_1} \tag{17}$$

### 8.2 Two-parameter case

Select a reference point $\boldsymbol{\mu}^\star = (\mu^\star, \beta^\star)^T$ and a unit direction
$\bar{\boldsymbol{\mu}} = (\bar{\mu}, \bar{\beta})^T$ with $|\bar{\boldsymbol{\mu}}| = 1$.

Write $\tilde{\boldsymbol{\mu}} = \boldsymbol{\mu}^\star + \tilde{\delta}\bar{\boldsymbol{\mu}}$ and
substitute into the fitting polynomial to obtain a scalar quadratic in $\tilde{\delta}$:

$$A_\delta \tilde{\delta}^2 + B_\delta \tilde{\delta} + C_\delta = 0 \tag{18}$$

where:

$$A_\delta = \bar{\boldsymbol{\mu}}^T \mathbf{a}_2 \bar{\boldsymbol{\mu}}$$

$$B_\delta = \mathbf{a}_1^T\bar{\boldsymbol{\mu}} + 2(\boldsymbol{\mu}^\star)^T\mathbf{a}_2\bar{\boldsymbol{\mu}}$$

$$C_\delta = a_0 + \mathbf{a}_1^T\boldsymbol{\mu}^\star + (\boldsymbol{\mu}^\star)^T\mathbf{a}_2\boldsymbol{\mu}^\star$$

and $\mathbf{a}_2$ is the $2 \times 2$ symmetric matrix:

$$\mathbf{a}_2 = \begin{bmatrix} a_2^{(\mu\mu)} & \tfrac{1}{2}a_2^{(\mu\beta)} \\ \tfrac{1}{2}a_2^{(\mu\beta)} & a_2^{(\beta\beta)} \end{bmatrix} \tag{19}$$

For the first-order fit, $\mathbf{a}_2 = \mathbf{0}$, reducing Eq. 18 to a linear equation
$\tilde{\delta} = -C_\delta/B_\delta$.

**Tracing the flutter boundary ($\tilde{r} = 0$):**

Fix $\bar{\boldsymbol{\mu}} = (1, 0)^T$ (direction of increasing $\mu$ at fixed $\beta$) and sweep
$\beta^\star$ from 0.25 to 1.75 in steps of 0.05. At each $\beta^\star$, set $\mu^\star = -0.5$ and
solve Eq. 18 for $\tilde{\delta}$; the flutter boundary point is $(\mu_c, \beta^\star) = \mu^\star + \tilde{\delta}\bar{\boldsymbol{\mu}}$.

**Tracing bifurcation diagram at fixed $\beta$:**

To obtain sections of the diagram at $\beta = 0.5, 1.0, 1.5$, fix $\bar{\boldsymbol{\mu}} = (1,0)^T$,
set $\beta^\star$ to the desired value, and vary $\mu^\star$ from $-0.5$ upward in small steps. At
each amplitude $\tilde{r}_k$, solve Eq. 18 to get $\tilde{\mu}_k(\beta^\star)$.

---

## 9. Expected Results and Verification Checks

### Check 1: Recovery rate accuracy (Step 5)

For each simulation $l$ and each time step $k$, compute:

$$\epsilon_\lambda^{(l,k)} = \left|\lambda_k^{\mathrm{FD}} - \lambda^{\mathrm{exact}}(r_k, \mu_l, \beta_l)\right|$$

where $\lambda^{\mathrm{exact}}$ is Eq. 3 or 4. The maximum error should satisfy:
- Supercritical system: $\max \epsilon_\lambda < 5 \times 10^{-3}$ for $\Delta t = 0.01$.
- Subcritical system: $\max \epsilon_\lambda < 7 \times 10^{-3}$; the $r^4$ term gives a steeper $\lambda(r)$ gradient at large $r$, amplifying the first-order FD truncation error.

### Check 2: Fitting exactness (Section 7)

The true recovery rate is exactly bi-linear in $(\mu, \beta)$ (Eqs. 3–4), so a first-order surface fit to the overdetermined system (9 curves, 3 unknown coefficients per amplitude) should reproduce each input curve to within the finite-difference noise floor:

$$\left|\hat{\lambda}(\tilde{r}_k, \mu_l, \beta_l) - \lambda_{kl}^{\mathrm{FD}}\right| < 10^{-3} \quad \forall\, k, l$$

where $\hat{\lambda}$ is the fitted surface and $\lambda_{kl}^{\mathrm{FD}}$ are the FD input curves. Machine-precision agreement ($10^{-10}$) is not achievable from an overdetermined LS problem fit to 9 noisy curves; the tolerance is governed by the FD truncation error $O(\Delta t^2) \approx 10^{-4}$, not the algebraic exactness of the system.

### Check 3: Flutter boundary

The predicted flutter boundary must satisfy $\mu_c \approx 0$ for all values of $\beta \in [0.25, 1.75]$.
Acceptable tolerance: $|\mu_c| < 10^{-3}$.

### Check 4: Bifurcation diagram at $\beta = 1$

Compare the predicted diagram against Eqs. 5–7 at $\beta = 1$ over $r \in [0, 1.4]$.

For the supercritical system ($\tilde{\mu}^{\mathrm{exact}} = \beta\tilde{r}^2$ at $\beta = \alpha = 1$):

$$\left|\tilde{\mu}_k^{\mathrm{pred}} - \tilde{\mu}^{\mathrm{exact}}(\tilde{r}_k)\right| < 2 \times 10^{-3}$$

The slack beyond the FD noise floor arises from the overdetermined LS fit; errors are largest near $r = r_{\max}$ where the transient is close to its initial condition. The bifurcation type should be identified as `:supercritical`.

For the subcritical system ($\tilde{\mu}^{\mathrm{exact}} = \tilde{r}^4 - \tilde{r}^2$ at $\beta = \gamma = 1$):

$$\left|\tilde{\mu}_k^{\mathrm{pred}} - \tilde{\mu}^{\mathrm{exact}}(\tilde{r}_k)\right| < 1.5 \times 10^{-2}$$

The wider tolerance reflects the $r^4$ gradient amplifying FD errors into predicted-$\mu$ errors at large $r$. The bifurcation type should be identified as `:subcritical`.

### Check 5: Two-parameter surface

The predicted bifurcation diagram surface in $(\mu, \beta, r)$-space should lie on the analytical
surfaces (Eqs. 5–7) to within $10^{-3}$ at points within the range of training data, and
extrapolate correctly to nearby untested $\beta$ values.

**Key cross-check:** Use samples 4–9 only (no samples at $\beta = 1$). The predicted diagram at
$\beta = 1$ should still match the analytical solution to within $10^{-2}$. This verifies that the
method correctly interpolates to untested parameter combinations.

### Check 6: Single- vs. two-parameter agreement

At $\beta = 1$, the one-parameter method (samples 1–3) and the two-parameter method (samples 1–9)
should give identical bifurcation diagrams, since the same data at $\beta = 1$ is used in both cases.

---

## 10. Illustrative Figures to Reproduce

The following figures from Riso JFS 2021 can be used as quantitative references:

- **Fig. 5 (supercritical, panels a–b and subcritical, panels c–d):** Transient responses and
  recovery rate curves for simulations 1–3. Recovery rate curves in the $r$--$\lambda$ plane should
  show monotonically decreasing $\lambda$ with $r$ (supercritical) and non-monotonic behaviour
  (subcritical).

- **Fig. 6 (recovery rates for simulations 4–9):** Similar curves for the two-parameter samples.

- **Fig. 7:** Bifurcation diagrams at $\beta = 1$ in the $\mu$--$r$ plane for both one- and
  two-parameter methods, compared with the analytical solution. All three should coincide.

- **Fig. 8:** Bifurcation diagram surfaces in the $\mu$--$\beta$--$r$ space (four panels:
  supercritical/subcritical, linear/quadratic fitting). Red markers show training samples; black
  markers show forecasted points; the surface is the analytical solution.

---

## 11. Implementation Notes

### ODE solver

For the supercritical system, the exact solution is available:

$$r(t) = r_0 \left[\frac{\alpha\mu}{\beta r_0^2 + (\alpha\mu - \beta r_0^2)e^{2\alpha\mu t}}\right]^{1/2}$$

This can be used to generate noise-free data without any ODE solver, which is useful for isolating
fitting errors from integration errors.

For general use, a simple RK4 integrator with step $\Delta t = 0.01$ is sufficient.

### Handling the non-oscillatory nature

Because these are scalar ODEs, $r(t)$ is the amplitude directly — there is no need for phase fixing,
bandpass filtering, or ERA. The full time series can be used in the finite-difference estimator
(Eq. 8) at every time step, not just at peaks.

### Initial transient

For these simple ODEs there is no initial transient from non-bifurcating modes. Begin using data
from $t = 0$.

### Edge cases

- At very small $r$ (near zero), the recovery rate estimate becomes noisy due to logarithm of small
  numbers. Discard data points with $r < r_{\min} = 0.01$.
- Near the end of the transient ($r \approx 0$), the signal may underflow. Terminate integration
  when $r < 10^{-6}$.
