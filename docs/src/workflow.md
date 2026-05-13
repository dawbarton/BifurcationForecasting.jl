# Workflow

This page walks through a complete bifurcation forecasting example using a
one-degree-of-freedom amplitude ODE with a known analytical solution.  The
system is the supercritical Hopf normal form

```
ṙ = r (μ − r²),    r ≥ 0
```

with control parameter `μ`.  The true flutter boundary is at `μ = 0`; the
stable LCO branch satisfies `r = √μ` for `μ > 0`.

The goal is to estimate the flutter boundary and LCO branch from free-decay
transients measured at three pre-flutter values `μ ∈ {−0.5, −0.4, −0.3}`.

## Step 1 — Generate synthetic transients

```julia
using BifurcationForecasting
using OrdinaryDiffEq

# Right-hand side of the amplitude ODE.
function rhs!(du, u, p, t)
    μ = p[1]
    du[1] = u[1] * (μ - u[1]^2)
end

# Simulate a free-decay transient at a given μ.
function simulate(μ; r0 = 1.5, dt = 0.05, t_max = 60.0)
    prob = ODEProblem(rhs!, [r0], (0.0, t_max), [μ])
    sol  = solve(prob, Tsit5(); saveat = dt, abstol = 1e-10, reltol = 1e-8)
    return Float64.(sol.t), Float64.(first.(sol.u))
end

μ_samples = [-0.5, -0.4, -0.3]
```

## Step 2 — Wrap data in `TransientSample` / `TransientDataset`

[`TransientSample`](@ref) holds a single transient.  The scalar field `q` must
already be a modal amplitude signal — for this 1-DOF example the amplitude
`r(t)` plays that role directly.

```julia
samples = TransientSample[]
for μ in μ_samples
    t, r = simulate(μ)
    push!(samples, TransientSample([μ], t, r))
end

dataset = TransientDataset(samples, ["μ"])
```

## Step 3 — Estimate recovery rates

For non-oscillatory systems (or after peak extraction), the
[`FiniteDifference`](@ref) estimator requires only consecutive amplitude pairs.
For oscillatory systems with few peaks, prefer [`EnvelopeODE`](@ref).

```julia
r_grid = range(0.05, 1.45; length = 50) |> collect

curves = RecoveryRateCurve[]
for s in dataset.samples
    # find_peaks / discard_initial_transient are not needed here because r(t) is
    # already a monotonically decaying amplitude envelope (not an oscillating signal).
    # For oscillatory data, call find_peaks → discard_initial_transient first.
    curve = estimate_recovery_rate(s.t, s.q, s.μ, FiniteDifference())
    push!(curves, curve)
end
```

!!! note "Oscillatory signals"
    For physically oscillatory systems (e.g. aeroelastic sections), first
    project onto the bifurcating mode with [`modal_project`](@ref) (after
    [`era_modal_basis`](@ref)) or [`bandpass_filter`](@ref), then call
    [`find_peaks`](@ref) and [`discard_initial_transient`](@ref) before
    estimating recovery rates.

## Step 4 — Fit the polynomial surface

[`fit_polynomial_surface`](@ref) fits a polynomial approximation `λ(μ, r)`
from the estimated curves.  With a single control parameter and `order = 1`
this is linear in `μ`.

```julia
surface = fit_polynomial_surface(curves, r_grid; order = 1)
```

## Step 5 — Check conditioning (optional)

```julia
κ = check_condition(curves, r_grid; order = 1)
println("Design-matrix condition number: ", κ[1])
```

A condition number above `1e6` suggests the parameter samples are poorly
spaced; consider adding more samples or enabling ridge regularisation via the
`ridge` keyword of `fit_polynomial_surface`.

## Step 6 — Extract the bifurcation diagram

[`extract_diagram`](@ref) finds the zero-recovery-rate locus `λ(μ, r) = 0` at
each amplitude grid point and classifies the bifurcation type.

```julia
diagram = extract_diagram(surface)

println("Flutter boundary: μ_c = ", diagram.flutter_boundary[1][1])
println("Bifurcation type: ",        diagram.bifurcation_type)

# Compare to the analytical LCO branch r = √μ.
for pt in diagram.points[1:5:end]
    μ_pred  = pt.μ[1]
    μ_exact = pt.r^2
    @printf "r = %.3f  μ_pred = %+.4f  μ_exact = %+.4f\n" pt.r μ_pred μ_exact
end
```

Expected output (approximate):

```
Flutter boundary: μ_c = 0.0002
Bifurcation type: supercritical
r = 0.050  μ_pred = +0.0025  μ_exact = +0.0025
r = 0.345  μ_pred = +0.1191  μ_exact = +0.1190
...
```

## High-level alternative: `forecast`

The [`forecast`](@ref) function wraps Steps 3–4 into a single call.  It
expects `TransientSample.q` to already be a modal amplitude signal and also
requires the flutter-mode angular frequency `ω_c` (needed to discard the
initial transient for oscillatory data).

```julia
# For oscillatory signals with known flutter frequency ω_c ≈ 6.28 rad/s:
# surface = forecast(dataset, ω_c; method = EnvelopeODE(), order = 1)

# For the non-oscillatory 1-DOF case we call the lower-level API directly
# (as shown above), since find_peaks / discard_initial_transient are bypassed.
```

## State velocity method

If the state equations `ẏ = f(μ, y)` and their Jacobian are known, recovery
rates can be computed *without* time integration using
[`forecast_state_velocity`](@ref):

```julia
# Vector-field and Jacobian for a 1-DOF system (y = [r]).
f(μ, y) = [y[1] * (μ[1] - y[1]^2)]
J(μ, y) = [μ[1] - 3*y[1]^2;;]   # 1×1 matrix

μ_vecs  = [[μ] for μ in μ_samples]
r_grid  = range(0.05, 1.45; length = 50) |> collect

surface_sv = forecast_state_velocity(f, J, μ_vecs, r_grid; n_states = 1)
diagram_sv = extract_diagram(surface_sv)
println("State-velocity flutter boundary: ", diagram_sv.flutter_boundary[1][1])
```

## Convergence check

When it is unclear whether a first-order polynomial surface is sufficient,
[`convergence_check`](@ref) compares first- and second-order fits:

```julia
result = convergence_check(curves, r_grid)
println("Converged (1st ≈ 2nd order): ", result.converged)
println("Relative difference in flutter speed: ", result.rel_diff)
```

## Latin Hypercube Sampling

For multi-parameter problems, [`lhs_samples`](@ref) generates well-distributed
parameter samples:

```julia
# 8 samples in a 2D parameter box: μ ∈ [−0.6, −0.2], β ∈ [0.8, 1.2].
μ_lhs = lhs_samples(8, [-0.6, 0.8], [-0.2, 1.2])
```

[`amplitude_grid`](@ref) provides square-root-spaced grids for better
resolution near the flutter boundary:

```julia
r_grid = amplitude_grid(50, 1.5; r_min = 0.01, spacing = :sqrt)
```
