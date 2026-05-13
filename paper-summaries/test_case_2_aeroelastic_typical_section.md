# Test Case 2: Two-DOF Aeroelastic Typical Section

*Based on Ghadami & Epureanu, J. Comput. Nonlinear Dyn. 11 (2016) 061009 (transient-data method);
Riso, Cesnik & Epureanu, AIAA Journal 60 (2022) 5401–5413 (state velocity method), Section 4.*

---

## 1. Purpose

This test case verifies the full bifurcation forecasting pipeline — including ERA modal filtering,
phase fixing, nonlinear optimisation for recovery rates, polynomial fitting, and bifurcation diagram
extraction — on a four-state oscillatory aeroelastic system. It also verifies the state velocity
extension. The flutter speed and bifurcation diagram shape are compared against eigenvalue analysis
and direct time-marching reference solutions.

---

## Note on Automated Test Implementation

The automated test suite (`test/test_aeroelastic_full.jl`) implements the full four-state typical
section described below. The tests exercise ERA modal filtering, biorthogonal modal projection,
`EnvelopeODE` recovery-rate estimation, polynomial surface fitting, bifurcation diagram extraction,
and the state velocity method against the actual system.

**Automated test parameters:**

| Quantity | Value |
|---|---|
| Parameter samples $\bar{U}$ (supercritical) | $\{0.825,\,0.830,\,0.835\}$ |
| ERA model order | 8 (4 conjugate pairs; 2 physical, 2 spurious — see ERA note below) |
| ERA Hankel dimensions | $r_{\mathrm{H}} = c_{\mathrm{H}} = 500$ |
| Recovery-rate estimator | `EnvelopeODE(4)` (order $p = 4$) |
| Surface order | 2 (quadratic in $\bar{U}$; exact with 3 samples) |
| Amplitude normalisation | peaks divided by first post-discard peak ($\tilde{r}_0 \equiv 1$) |
| Initial transient discard | $n_{\mathrm{periods}} = 5$ periods of $\omega_c$ |
| Simulation: step and duration | $\Delta\bar{t} = 0.5$, $\bar{t}_{\max} = 1000$ |
| State velocity $r$-grid | 40 points in $[10^{-4},\,0.3]$ (physical units, displacement along $\mathrm{Re}(\boldsymbol{\phi}_c)$) |

**Automated test checks (7 total):**

| Check | Quantity | Tolerance |
|---|---|---|
| 1 | Flutter speed from eigenvalue sweep: $\bar{U}_F = 0.832$ | 0.001 |
| 2 | Reduced modal matrix $\tilde{\mathbf{A}}_c$ recovers $g_c$ and $\omega_c$ | $10^{-8}$ (matrix entries) |
| 3 | $\lambda(r \to 0) \approx g_c(\bar{U})$ from `EnvelopeODE` | 0.05 |
| 4 | Flutter boundary from polynomial surface within $|{}\cdot{} - 0.832|$ | 0.005 |
| 5 | Bifurcation type is `:supercritical` | exact |
| 6 | All LCO branch points above $\bar{U}_F$ | $\bar{U}_F - 0.005$ |
| 7 | State velocity surface: $g_c$ recovery and flutter boundary | 0.001 / 0.005 |

**Note on Check 7:** The transient-data and state velocity surfaces use incompatible amplitude
normalisations — ERA-projected peak amplitudes versus physical displacement along $\mathrm{Re}(\boldsymbol{\phi}_c)$
— so a direct pointwise $\lambda(r)$ comparison requires a known scale factor and is not implemented.
Instead Check 7 verifies that the state velocity surface recovers $g_c(\bar{U})$ at small $r$ (the
limit where both methods must agree) and predicts the flutter boundary within 0.005.

---

---

## 2. System Description

A two-degree-of-freedom typical section (pitch and plunge) with quasi-steady incompressible
aerodynamics and polynomial structural nonlinearities in pitch. The nondimensionalised state vector
is $\mathbf{y} = \{\bar{h}, \alpha, \dot{\bar{h}}, \dot{\alpha}\}^T$ where $\bar{h} = h/b$ is the
normalised plunge displacement and $\alpha$ is the pitch angle in radians.

---

## 3. Governing Equations

### 3.1 Dimensional equations of motion

$$m\ddot{h} + S_\alpha\ddot{\alpha} + K_h^{(1)}h = -2\pi\rho b U^2\!\left(\alpha + \frac{\dot{h}}{U}\right) \tag{1}$$

$$S_\alpha\ddot{h} + I_\alpha\ddot{\alpha} + K_\alpha^{(1)}\alpha + K_\alpha^{(3)}\alpha^3 + K_\alpha^{(5)}\alpha^5 = 2\pi\rho b U^2\!\left(\alpha + \frac{\dot{h}}{U}\right)e \tag{2}$$

### 3.2 Nondimensionalisation

Define nondimensional parameters:

$$\bar{m} = \frac{m}{\rho\pi b^2}, \quad x_\alpha = \frac{S_\alpha}{mb}, \quad \Omega^2 = \frac{\omega_h^2}{\omega_\alpha^2}, \quad r_\alpha^2 = \frac{I_\alpha}{mb^2}, \quad \bar{e} = \frac{e}{b}$$

$$\bar{U} = \frac{U}{b\omega_\alpha}, \quad \bar{t} = t\omega_\alpha, \quad \bar{h} = \frac{h}{b}, \quad \kappa_\alpha^{(3)} = \frac{K_\alpha^{(3)}}{K_\alpha^{(1)}}, \quad \kappa_\alpha^{(5)} = \frac{K_\alpha^{(5)}}{K_\alpha^{(1)}}$$

where $\omega_h^2 = K_h^{(1)}/m$ and $\omega_\alpha^2 = K_\alpha^{(1)}/I_\alpha$.

### 3.3 Nondimensional state-space form

$$\dot{\mathbf{y}} = \mathbf{f}(\bar{U}, \mathbf{y}) = \mathbf{A}(\bar{U})\mathbf{y} + \mathbf{f}_{nl}(\mathbf{y}) \tag{3}$$

with state vector $\mathbf{y} = \{\bar{h}, \alpha, \dot{\bar{h}}, \dot{\alpha}\}^T$.

### 3.4 Linear state matrix

$$\mathbf{A}(\bar{U}) = \begin{bmatrix} \mathbf{0}_{2\times2} & \mathbf{I}_{2\times2} \\ -\mathbf{M}_S^{-1}[\mathbf{K}_S + \mathbf{K}_A(\bar{U})] & -\mathbf{M}_S^{-1}\mathbf{D}_A(\bar{U}) \end{bmatrix} \tag{4}$$

**Structural mass matrix:**

$$\mathbf{M}_S = \begin{bmatrix} 1 & x_\alpha \\ x_\alpha & r_\alpha^2 \end{bmatrix} \tag{5}$$

**Structural stiffness matrix:**

$$\mathbf{K}_S = \begin{bmatrix} \Omega^2 & 0 \\ 0 & r_\alpha^2 \end{bmatrix} \tag{6}$$

**Aerodynamic damping matrix:**

$$\mathbf{D}_A(\bar{U}) = \frac{2}{\bar{m}}\bar{U}\begin{bmatrix} 1 & 0 \\ -\bar{e} & 0 \end{bmatrix} \tag{7}$$

**Aerodynamic stiffness matrix:**

$$\mathbf{K}_A(\bar{U}) = \frac{2}{\bar{m}}\bar{U}^2\begin{bmatrix} 0 & 1 \\ 0 & -\bar{e} \end{bmatrix} \tag{8}$$

### 3.5 Nonlinear forcing

$$\mathbf{f}_{nl}(\mathbf{y}) = \frac{\kappa_\alpha^{(3)}\alpha^3 + \kappa_\alpha^{(5)}\alpha^5}{1 - x_\alpha^2/r_\alpha^2}\{0,\, 0,\, x_\alpha,\, -1\}^T \tag{9}$$

This vector is derived from the inverse of $\mathbf{M}_S$ applied to the nonlinear pitch restoring
force. The denominator $1 - x_\alpha^2/r_\alpha^2$ arises from the $2\times 2$ determinant of
$\mathbf{M}_S$.

**Explicit expansion of Eq. 9.** Let $\Delta = 1 - x_\alpha^2/r_\alpha^2$ and
$F_{nl} = \kappa_\alpha^{(3)}\alpha^3 + \kappa_\alpha^{(5)}\alpha^5$. Then:

$$f_{nl,1} = 0, \quad f_{nl,2} = 0, \quad f_{nl,3} = \frac{x_\alpha F_{nl}}{\Delta}, \quad f_{nl,4} = \frac{-F_{nl}}{\Delta} \tag{10}$$

---

## 4. System Parameters

### 4.1 Linear parameters

| Parameter | Value |
|---|---|
| Mass ratio $\bar{m}$ | 10 |
| Nondimensional static unbalance $x_\alpha$ | 0.2 |
| Nondimensional radius of gyration $r_\alpha$ | 0.3 |
| Nondimensional elastic axis offset $\bar{e}$ | 0.2 |
| Frequency ratio $\Omega = \omega_h/\omega_\alpha$ | 0.5 |

### 4.2 Nonlinear parameters

| Case | $\kappa_\alpha^{(3)}$ | $\kappa_\alpha^{(5)}$ | Bifurcation type |
|---|---|---|---|
| Supercritical | 1.5 | 0 | Supercritical Hopf |
| Subcritical | −1.5 | 50 | Subcritical Hopf |

---

## 5. Linear Analysis

### 5.1 Flutter speed

Compute eigenvalues $\sigma_l(\bar{U}) = g_l(\bar{U}) \pm j\omega_l(\bar{U})$ of $\mathbf{A}(\bar{U})$
for a range of $\bar{U}$ values (e.g.\ $\bar{U} \in [0, 1.4]$ in steps of 0.01).

The **flutter speed** is the value $\bar{U}_F$ where the real part of one eigenvalue first crosses
zero:

$$g_c(\bar{U}_F) = 0 \tag{11}$$

**Expected result:** $\bar{U}_F = 0.832$.

### 5.2 Bifurcating mode identification

At $\bar{U}$ slightly below $\bar{U}_F$, the bifurcating mode is the one with the most positive (least
negative) real part. Extract its:
- Eigenvalue: $\sigma_c = g_c + j\omega_c$
- Right eigenvector: $\boldsymbol{\phi}_c \in \mathbb{C}^4$
- Left eigenvector: $\boldsymbol{\psi}_c \in \mathbb{C}^4$

The left eigenvector satisfies $\boldsymbol{\psi}_c^T \mathbf{A} = \sigma_c \boldsymbol{\psi}_c^T$
(equivalently, $\mathbf{A}^T\boldsymbol{\psi}_c^* = \sigma_c^*\boldsymbol{\psi}_c^*$).

Form the $4 \times 2$ matrices:

$$\boldsymbol{\Phi} = \begin{bmatrix} \mathrm{Re}(\boldsymbol{\phi}_c) & \mathrm{Im}(\boldsymbol{\phi}_c) \end{bmatrix}, \quad \boldsymbol{\Psi} = \begin{bmatrix} \mathrm{Re}(\boldsymbol{\psi}_c) & \mathrm{Im}(\boldsymbol{\psi}_c) \end{bmatrix} \tag{12}$$

**Verification:** The reduced $2\times 2$ matrix should recover the flutter mode eigenvalue:

$$\tilde{\mathbf{A}}_c = (\boldsymbol{\Psi}^T\boldsymbol{\Phi})^{-1}\boldsymbol{\Psi}^T\mathbf{A}\boldsymbol{\Phi} = \begin{bmatrix} g_c & \omega_c \\ -\omega_c & g_c \end{bmatrix} \tag{13}$$

Check that the off-diagonal entries of $\tilde{\mathbf{A}}_c$ satisfy
$|\tilde{A}_{12} + \tilde{A}_{21}| < 10^{-10}$ and $|\tilde{A}_{12} - \omega_c| < 10^{-8}$.

---

## 6. Transient-Data Method

### 6.1 Perturbations and pre-flutter simulations

Apply initial conditions along the real part of the bifurcating mode eigenvector:

$$\mathbf{y}(0) = \boldsymbol{\Phi}\begin{pmatrix}r_0 \\ 0\end{pmatrix} = r_0\,\mathrm{Re}(\boldsymbol{\phi}_c) \tag{14}$$

with $r_0 = 15°$ (convert to radians for $\alpha$; scale $\bar{h}$ component proportionally).

Integrate Eq. 3 forward in time using RK4 with nondimensional time step $\Delta\bar{t} = 0.5$
(corresponding to roughly $\omega_c \Delta\bar{t}/(2\pi) \approx 0.08$ steps per cycle) until
$t_{\max} = 1000$ or until $|\alpha| < 10^{-4}$ rad.

**Parameter samples for transient-data method:**

| Simulation | $\bar{U}$ | Use |
|---|---|---|
| 1 | 0.825 | Pre-flutter (supercritical) |
| 2 | 0.830 | Pre-flutter (supercritical) |
| 3 | 0.835 | Pre-flutter (supercritical) |

For the subcritical case, use $\bar{U} = 0.820, 0.825, 0.830$.

### 6.2 Modal projection

Project each transient response onto the bifurcating mode using the $\boldsymbol{\Phi}$,
$\boldsymbol{\Psi}$ matrices computed in Section 5.2 (evaluated at $\bar{U}$ closest to flutter):

The projected modal coordinate is obtained by applying the biorthogonal projection:

$$\mathbf{q}(t) = (\boldsymbol{\Psi}^T\boldsymbol{\Phi})^{-1}\boldsymbol{\Psi}^T\Delta\mathbf{y}(t) \tag{15}$$

where $\Delta\mathbf{y}(t) = \mathbf{y}(t) - \mathbf{y}_e$ is the perturbation from equilibrium
(here $\mathbf{y}_e = \mathbf{0}$ since there is no mean flow angle). The projected scalar signal
(for phase $\theta = 0$) is:

$$q(t) = q_1(t) = \text{first component of } \mathbf{q}(t) \tag{16}$$

This signal oscillates at frequency $\approx \omega_c$ and decays at rate $\approx g_c$.

**Alternative (simpler):** Since the pitch angle $\alpha$ is the dominant output, use it directly
after bandpass filtering. The filter should have:
- Centre frequency: $\omega_c / (2\pi)$ (in the nondimensional frequency units of $\bar{t}$)
- Bandwidth: $\pm 20\%$ of centre frequency

For this system, $\omega_c \approx 0.83$ rad per nondimensional time unit (see Section 5), so the
filter passes frequencies in the range $[0.66, 1.00]$ nondimensional rad per time unit.

### 6.3 Phase fixing

Extract local maxima of $|q(t)|$ (or $|\alpha(t)|$ if using the bandpass approach):

1. Find all time indices $k$ where $q(t_k) > q(t_{k-1})$ and $q(t_k) > q(t_{k+1})$.
2. Retain only positive maxima (phase $\theta = 0^+$): $r_k = q(t_k) > 0$.
3. Discard the first $T_{\mathrm{skip}} = 5$ oscillation periods ($\approx 5 \times 2\pi/\omega_c$
   nondimensional time units) to allow non-bifurcating modes to decay.

### 6.4 Recovery rate estimation

Because the oscillation frequency is low relative to the decay rate at parameter values far from
flutter, use the **nonlinear optimisation approach** (always preferred over finite differences for
oscillatory systems).

At each parameter sample $\bar{U}_i$, fit the envelope ODE:

$$\dot{r} = r\left(\lambda_0 + \lambda_1 r + \lambda_2 r^2 + \lambda_3 r^3 + \lambda_4 r^4\right) \tag{17}$$

to the observed peaks $(t_k, r_k)$.

**Optimisation procedure:**

1. **Parameterise:** Unknown vector $\boldsymbol{\lambda} = [\lambda_0, \lambda_1, \lambda_2, \lambda_3, \lambda_4]^T$.
2. **Forward integration:** Given $\boldsymbol{\lambda}$ and initial condition $r(t_0) = r_0^{\mathrm{obs}}$,
   integrate Eq. 17 using RK4 with step $\Delta t_{\mathrm{env}} = 0.1$ (coarser than full simulation
   since the envelope varies slowly).
3. **Residual:** At each peak time $t_k$, evaluate $r^{\mathrm{model}}(t_k; \boldsymbol{\lambda})$ and
   compute the residual $\varepsilon_k = r^{\mathrm{model}}(t_k) - r_k^{\mathrm{obs}}$.
4. **Minimise:** Use Levenberg–Marquardt (or any nonlinear least squares solver) to minimise
   $\sum_k \varepsilon_k^2$.
5. **Initialisation:** Start from $\boldsymbol{\lambda}^{(0)} = [g_c, 0, 0, 0, 0]^T$ where $g_c$ is
   the linear damping of the bifurcating mode at $\bar{U}_i$.

Once $\boldsymbol{\lambda}$ is found for a given $\bar{U}_i$, the recovery rate at any amplitude is:

$$\lambda(\bar{U}_i, r) = \lambda_0 + \lambda_1 r + \cdots + \lambda_4 r^4 \tag{18}$$

**Order selection:** Run with $p = 2, 3, 4$ and compare predicted flutter speeds; $p = 4$ should
give a stable result. Use $p = 4$ as the default.

### 6.5 Amplitude grid

Select $N_r = 40$ amplitudes uniformly spaced in $[0.1°, 12°]$ (in degrees for pitch angle, or
equivalently in the nondimensional modal amplitude $r$).

For each amplitude $\tilde{r}_k$ and each simulation $i$, evaluate $\lambda(\bar{U}_i, \tilde{r}_k)$
from the fitted polynomial (Eq. 18).

### 6.6 Polynomial fitting in $\bar{U}$

At each amplitude $\tilde{r}_k$, fit a second-order polynomial in $\bar{U}$ (Eq. 9 of the primary
algorithm document, with $\mu \to \bar{U}$) using the three data points:

$$\mathbf{A} = \begin{bmatrix} 1 & \bar{U}_1 & \bar{U}_1^2 \\ 1 & \bar{U}_2 & \bar{U}_2^2 \\ 1 & \bar{U}_3 & \bar{U}_3^2 \end{bmatrix}, \quad \boldsymbol{\lambda}_k = \begin{bmatrix} \lambda(\bar{U}_1, \tilde{r}_k) \\ \lambda(\bar{U}_2, \tilde{r}_k) \\ \lambda(\bar{U}_3, \tilde{r}_k) \end{bmatrix} \tag{19}$$

Solve: $\mathbf{c}_k = \mathbf{A}^{-1}\boldsymbol{\lambda}_k$ (system is square; use LU factorisation).

### 6.7 Bifurcation diagram extraction

At each amplitude $\tilde{r}_k$, find $\tilde{\bar{U}}_k$ such that:

$$a_0(\tilde{r}_k) + a_1(\tilde{r}_k)\tilde{\bar{U}}_k + a_2(\tilde{r}_k)\tilde{\bar{U}}_k^2 = 0 \tag{20}$$

Use the quadratic formula; select the larger root (the one above the pre-flutter samples).

**Flutter speed estimate:** Evaluate at $\tilde{r} \to 0$ (use the smallest available amplitude
$\tilde{r}_1$). The predicted flutter speed should satisfy $|\tilde{\bar{U}}_1 - 0.832| < 0.005$.

---

## 7. State Velocity Method

This section describes the model-based variant that bypasses transient simulation. It requires the
state-space model (Section 3) directly.

### 7.1 Setup

Choose $N_p = 3$ pre-flutter flow speed values, e.g.\ $\bar{U}_i \in \{0.825, 0.830\}$ (two values
suffice for a linear fit in $\bar{U}$; three for a quadratic fit).

Choose an amplitude grid $r_k$ spanning $[0.001, 0.20]$ in nondimensional units (corresponding to
pitch angles $\approx 0.06°$ to $\approx 11.5°$ for the pitch component of the eigenvector).

### 7.2 Per-sample procedure

For each $\bar{U}_i$:

1. **Equilibrium:** The equilibrium is $\mathbf{y}_e = \mathbf{0}$ (zero angle of attack, no gravity).

2. **Jacobian:** Compute $\mathbf{A}(\bar{U}_i)$ from Eqs. 4–8.

3. **Eigenvectors:** Find right eigenvector $\boldsymbol{\phi}_{ic}$ and left eigenvector
   $\boldsymbol{\psi}_{ic}$ of the bifurcating mode (most positive real part). Form $\boldsymbol{\Phi}_{ic}$
   and $\boldsymbol{\Psi}_{ic}$ as in Eq. 12.

4. **Reduced matrix:** Compute $\tilde{\mathbf{A}}_{ic}$ via Eq. 13 to obtain $g_{ic}$ and $\omega_{ic}$.

5. **Reduced eigenvectors:** The reduced right eigenvector is $\tilde{\boldsymbol{\phi}}^R_{ic} = \{1, 0\}^T$
   and reduced left eigenvector is $\tilde{\boldsymbol{\psi}}^R_{ic} = \{1, 0\}^T$ (these are fixed
   by the choice of $\boldsymbol{\Phi}$, $\boldsymbol{\Psi}$ normalisation via Eq. 13).

6. **Perturbed state:** For each amplitude $r_k$, construct:
   $$\mathbf{y}_{0}(r_k) = \boldsymbol{\Phi}_{ic}\begin{pmatrix}r_k \\ 0\end{pmatrix} = r_k\,\mathrm{Re}(\boldsymbol{\phi}_{ic}) \tag{21}$$

7. **State velocity:** Evaluate the full nonlinear vector field at the perturbed state:
   $$\mathbf{f}(\bar{U}_i, \mathbf{y}_0(r_k)) = \mathbf{A}(\bar{U}_i)\mathbf{y}_0(r_k) + \mathbf{f}_{nl}(\mathbf{y}_0(r_k)) \tag{22}$$

8. **Nonlinear part:** Extract:
   $$\mathbf{f}_{nl}(\mathbf{y}_0) = \mathbf{f}(\bar{U}_i, \mathbf{y}_0) - \mathbf{A}(\bar{U}_i)\mathbf{y}_0 \tag{23}$$

9. **Recovery rate:** Project via Eq. 20 of Riso AIAA 2022:
   $$\lambda_{ik} = g_{ic} + \frac{(\tilde{\boldsymbol{\psi}}^R_{ic})^T(\boldsymbol{\Psi}_{ic}^T\boldsymbol{\Phi}_{ic})^{-1}\boldsymbol{\Psi}_{ic}^T\,\mathbf{f}_{nl}(\mathbf{y}_0(r_k))}{r_k} \tag{24}$$

   Expanded with $\tilde{\boldsymbol{\psi}}^R_{ic} = \{1,0\}^T$, Eq. 24 simplifies to:

   $$\lambda_{ik} = g_{ic} + \frac{\left[(\boldsymbol{\Psi}_{ic}^T\boldsymbol{\Phi}_{ic})^{-1}\boldsymbol{\Psi}_{ic}^T\,\mathbf{f}_{nl}(\mathbf{y}_0(r_k))\right]_1}{r_k} \tag{25}$$

   where $[\cdot]_1$ denotes the first component of the vector.

### 7.3 Fitting and diagram extraction

Identical to Steps 6.6–6.7 of the transient-data method, using the recovery rates $\lambda_{ik}$
from Eq. 25 instead of those from the nonlinear optimisation.

---

## 8. Reference Solutions

### 8.1 Direct time-marching (bifurcation diagram reference)

To generate reference bifurcation diagrams for comparison:

1. Choose post-flutter flow speeds $\bar{U} \in [0.832, 0.86]$ in steps of 0.005.
2. At each $\bar{U}$, apply an initial condition along $\mathrm{Re}(\boldsymbol{\phi}_c)$ with
   $r_0 = 15°$.
3. Integrate Eq. 3 forward until oscillations converge to a steady LCO (typically $\bar{t} = 2000$).
4. Record the LCO amplitude as the maximum of $|\alpha(t)|$ over the final 10 oscillation cycles.

This gives the reference $(\bar{U}, \tilde{r})$ diagram.

### 8.2 Expected flutter speed

$$\bar{U}_F = 0.832$$

### 8.3 Expected bifurcation types

- **Supercritical ($\kappa_\alpha^{(3)} = 1.5$):** Stable LCO grows from zero amplitude at
  $\bar{U}_F$. LCO amplitude at $\bar{U} = 0.85$ is approximately $10°$–$12°$ pitch.
- **Subcritical ($\kappa_\alpha^{(3)} = -1.5$, $\kappa_\alpha^{(5)} = 50$):** Unstable LCO branch
  below $\bar{U}_F$; stable large-amplitude LCO above. Bi-stability region is narrow
  (approximately $0.825 < \bar{U} < 0.832$). LCO amplitude at $\bar{U} = 0.85$ is approximately
  $12°$–$15°$ pitch.

---

## 9. Verification Checks

### Check 1: Flutter speed from eigenvalue analysis

Sweep $\bar{U}$ from 0 to 1.4 in steps of 0.001. Plot the real parts $g_l(\bar{U})$ of all four
eigenvalues. The first crossing of $g = 0$ should occur at $\bar{U}_F = 0.832 \pm 0.001$.

### Check 2: Modal projection validity

At $\bar{U} = 0.830$, compute $\tilde{\mathbf{A}}_c$ via Eq. 13. Verify:

- $\tilde{A}_{11} \approx \tilde{A}_{22} \approx g_c \approx -0.004$ (small negative number)
- $\tilde{A}_{12} \approx -\tilde{A}_{21} \approx \omega_c \approx 0.83$ (nondimensional flutter frequency)
- Off-diagonal symmetry: $|\tilde{A}_{12} + \tilde{A}_{21}| < 10^{-10}$

**Note on flutter frequency:** $\omega_c \approx 0.83$ rad per nondimensional time unit at
$\bar{U} = 0.830$. This is *not* the structural plunge frequency $\Omega = 0.5$; it is the flutter
frequency of the *coupled* aeroelastic mode, which lies between the two structural natural
frequencies (the uncoupled pitch-plunge system has natural frequencies $\approx 0.47$ and $1.42$
rad per nondimensional time unit at $\bar{U} = 0$).

### Check 3: Recovery rates from transient data

For the supercritical case at $\bar{U} = 0.825$, the recovery rate should:
- Satisfy $\lambda(r \to 0) \approx g_c(\bar{U} = 0.825) < 0$ (consistent with linear damping).
- Decrease monotonically with $r$ (more negative at larger amplitudes).
- Pass through approximately $\lambda \approx g_c$ at $r = 0$ and become more negative as $r$ increases.

For the subcritical case at $\bar{U} = 0.825$, the recovery rate should:
- Start at $\lambda \approx g_c < 0$ at $r = 0$.
- Increase with $r$ initially (less negative), reaching a local maximum near the unstable LCO amplitude.
- Decrease again at larger $r$.

### Check 4: Flutter speed from bifurcation forecasting

Using the transient-data method with 3 samples (supercritical case):

$$\left|\tilde{\bar{U}}_F^{\mathrm{pred}} - 0.832\right| < 0.005$$

Using the state velocity method with 2 samples:

$$\left|\tilde{\bar{U}}_F^{\mathrm{SV}} - 0.832\right| < 0.005$$

### Check 5: Bifurcation type identification

The sign of $\partial\lambda/\partial r|_{\lambda=0,\,r=0^+}$ distinguishes bifurcation type:
- Supercritical: $\partial\lambda/\partial r < 0$ at $r \to 0^+$ (recovery rate decreases with amplitude)
- Subcritical: $\partial\lambda/\partial r > 0$ at $r \to 0^+$ (recovery rate increases initially)

### Check 6: Bifurcation diagram shape

Compare the forecasted diagram with the direct time-marching reference (Section 8.1):

For the supercritical case, the forecasted diagram should satisfy:
- Flutter speed within 0.5% of $\bar{U}_F = 0.832$.
- LCO amplitude at $\bar{U} = 0.84$: within 20% of the time-marching value.
- Monotonically increasing amplitude with $\bar{U}$.

For the subcritical case, the forecasted stable branch (large amplitude) should broadly agree with
time-marching. Quantitative accuracy is lower for the unstable branch.

### Check 7: State velocity vs. transient-data comparison (Riso AIAA 2022 Fig. 10–11)

The state velocity method recovery rates should:
- Match the transient-data recovery rates at $r \to 0$ (both must converge to $g_c$).
- Capture the correct qualitative trend (monotonically decreasing for supercritical; non-monotonic
  for subcritical).
- Deviate from the transient-data rates at large amplitudes (this is expected).

The state velocity bifurcation diagrams (Fig. 11 of Riso 2022) should:
- Capture the correct bifurcation type in both cases.
- Give a flutter speed within 0.5% of 0.832.
- Show increasing deviation from time-marching at larger post-flutter amplitudes.

**Note on amplitude normalisation:** The transient-data and state velocity methods use
incompatible amplitude axes. The transient-data method normalises peak amplitudes by the first
post-discard peak of the ERA-projected signal $q(t)$ (an arbitrary scalar multiple of the physical
amplitude, set by the ERA eigenvector normalisation). The state velocity method uses physical
displacement along $\mathrm{Re}(\boldsymbol{\phi}_c)$, where $\boldsymbol{\phi}_c$ is the
unit-norm eigenvector from `eigen`. A direct pointwise comparison of $\lambda(r)$ curves therefore
requires multiplying the state velocity amplitude axis by the ERA projection scale factor (the value
of the first post-discard peak of $q(t)$ in physical units). Without this scale factor the two
$r$-axes are incommensurable. The limit $r \to 0$ is unaffected (both converge to $g_c$).

---

## 10. Illustrative Figures to Reproduce

The following published figures provide quantitative reference curves:

- **Riso AIAA 2022, Fig. 8:** Eigenvalue root locus and damping/frequency evolution vs. $\bar{U}$.
  The flutter crossing at $\bar{U} = 0.832$ should be clearly visible.

- **Riso AIAA 2022, Fig. 9:** Pre- and post-flutter transient responses at several $\bar{U}$ values
  for both cases. The envelope decay rate should visibly slow as $\bar{U} \to 0.832$.

- **Riso AIAA 2022, Fig. 10:** Recovery rate vs. pitch amplitude $\alpha$ (degrees) at two pre-flutter
  speeds for each case. The state velocity method (dots) should track the time-marching reference
  (solid line) at small amplitudes and diverge at large amplitudes.

- **Riso AIAA 2022, Fig. 11:** Bifurcation diagrams. Three curves per panel (time marching,
  bifurcation forecasting from transients, state velocity method). Supercritical panel: all three
  should closely agree. Subcritical panel: forecasting and time-marching agree; state velocity
  captures the type but diverges quantitatively.

---

## 11. Implementation Notes

### ODE solver recommendation

Use a fixed-step RK4 integrator with $\Delta\bar{t} = 0.5$ for generating transient responses. The
flutter frequency is $\omega_c \approx 0.83$ rad per nondimensional time unit (see Check 2 note),
so this gives approximately $2\pi/(0.83 \times 0.5) \approx 15$ steps per oscillation cycle —
adequate resolution.

For the envelope ODE (Eq. 17) in the nonlinear optimisation step, use $\Delta t_{\mathrm{env}} = 2$
(3–4 steps per oscillation cycle is sufficient for the slowly varying envelope).

### ERA implementation

If implementing ERA from scratch:

1. Assemble the multi-channel response matrix $Y_k \in \mathbb{R}^{n_{\mathrm{out}} \times n_{\mathrm{IC}}}$
   where $n_{\mathrm{out}} = 4$ (all state components) and $n_{\mathrm{IC}} = 1$ (single initial
   condition perturbation).
2. Choose Hankel dimensions $r = v = 500$ (as used by Ghadami & Epureanu 2016 for this system).
3. Perform SVD; truncate to the $N = 8$ largest singular values (4 modes, each with conjugate pair).
4. The bifurcating mode has eigenvalue with real part closest to zero (most positive real part
   among modes with positive imaginary part).

**Spurious modes with model_order = 8:** A 4-state system has only 2 physical conjugate pairs
(4 eigenvalues), but truncating the ERA SVD at 8 singular values introduces 2 additional spurious
modes. These arise from fitting noise or numerical rank in the Hankel matrix and have no physical
meaning. Selecting the mode with the most positive real part is robust to this: spurious ERA modes
typically have strongly negative real parts (they represent fast-decaying numerical artefacts) and
do not compete with the bifurcating mode near flutter.

**Simpler alternative for this system:** Since the system has only 4 states, skip ERA entirely and
compute eigenvectors directly from $\mathbf{A}(\bar{U})$ at the pre-flutter $\bar{U}$ closest to
flutter. This is valid for the state velocity method and for the transient-data method when a model
is available.

### Normalisation in optimisation

Scale the pitch angle to degrees before fitting the envelope ODE. Use $r$ in degrees throughout the
recovery rate estimation; convert back to radians only when evaluating $\mathbf{f}_{nl}$.

### Non-dimensionalisation check

To verify the nondimensionalisation, the reduced flutter frequency for this system is
$k_F = \omega_c/\bar{U}_F \approx 0.83/0.832 \approx 1.00$.
This is a physically reasonable value (reduced frequencies of order unity are typical for pitch-plunge
flutter in the incompressible regime).

### Constructing the nonlinear forcing vector (Eq. 9) explicitly

Given $\mathbf{y}_0 = \{y_1, y_2, y_3, y_4\}^T = \{\bar{h}, \alpha, \dot{\bar{h}}, \dot{\alpha}\}^T$:

$$\Delta = 1 - x_\alpha^2/r_\alpha^2 = 1 - 0.04/0.09 = 5/9$$

$$F_{nl}(\alpha) = \kappa_\alpha^{(3)}\alpha^3 + \kappa_\alpha^{(5)}\alpha^5$$

$$f_{nl,3} = \frac{x_\alpha}{\Delta}F_{nl}(\alpha) = \frac{0.2}{5/9}F_{nl} = 0.36\,F_{nl}$$

$$f_{nl,4} = \frac{-1}{\Delta}F_{nl}(\alpha) = -\frac{9}{5}F_{nl} = -1.8\,F_{nl}$$

$$f_{nl,1} = f_{nl,2} = 0$$
