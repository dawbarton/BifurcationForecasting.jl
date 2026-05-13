# BifurcationForecasting.jl

> **AI-generated code.** This package was written by Claude (Anthropic) under
> the direction of David A. W. Barton.  All scientific content has been
> reviewed by the author, but the implementation should be treated with
> appropriate scepticism and validated independently before use in research.

A Julia package for data-driven bifurcation forecasting based on the
*critical-slowing-down* (CSD) phenomenon.

Given a small set of free-decay transients measured in the pre-flutter regime,
the package reconstructs the full bifurcation diagram — flutter boundary and
limit-cycle oscillation (LCO) branch — **without crossing into the unstable
regime and without requiring a system model**.

A model-based variant (state velocity method) is also provided for cases where
the nonlinear state equations are available.

## Method overview

As a control parameter (e.g. airspeed) approaches a Hopf bifurcation, the
dominant eigenvalue approaches the imaginary axis and free-decay transients
slow down.  The *recovery rate*

```
λ(μ, r) = ṙ / r
```

quantifies this slow-down as a function of control parameter `μ` and oscillation
amplitude `r`.  The zero-level set `λ = 0` defines the bifurcation diagram:

- **Flutter boundary**: the locus `r = 0` where linear instability occurs.
- **LCO branch**: finite-amplitude oscillations at `r > 0`.

Recovery-rate curves are estimated from free-decay transients, interpolated onto
a common amplitude grid, and fitted with a polynomial in parameter space.
The diagram is then extracted by solving `λ(μ, r) = 0` analytically at each
amplitude.

## Primary references

- Ghadami & Epureanu (2016) *J. Comput. Nonlinear Dyn.* **11**, 061009.
  <https://doi.org/10.1115/1.4033920>
- Riso, Cesnik & Epureanu (2021) *J. Fluids Struct.* **101**, 103201.
  <https://doi.org/10.1016/j.jfluidstructs.2020.103201>
- Riso, Cesnik & Epureanu (2022) *AIAA J.* **60**, 5401–5413.
  <https://doi.org/10.2514/1.J061860>

## Installation

The package is not registered.  From the Julia REPL:

```julia
using Pkg
Pkg.add(url = "https://github.com/dawbarton/BifurcationForecasting.jl")
```

## Quick start

```julia
using BifurcationForecasting

# Wrap pre-processed transient data (one sample per parameter point).
# `q` must already be a modal amplitude signal: modal-projected or
# bandpass-filtered to isolate the bifurcating mode.
samples = [TransientSample([μ], t_vec, q_vec) for (μ, t_vec, q_vec) in data]
dataset = TransientDataset(samples, ["airspeed"])

# Run the full pipeline.  ω_c is the flutter-mode angular frequency (rad/s),
# used to determine how much of the initial transient to discard.
surface = forecast(dataset, ω_c; method = EnvelopeODE(), order = 1)

# Extract the bifurcation diagram.
diagram = extract_diagram(surface)
println("Flutter boundary: ", diagram.flutter_boundary)
println("Bifurcation type: ", diagram.bifurcation_type)   # :supercritical or :subcritical
```

## Typical workflow

```
Raw multi-channel data
        │
        ▼ era_modal_basis          (extract bifurcating mode via ERA)
        │ modal_project            (project data onto scalar modal coordinate)
        │   — or —
        │ bandpass_filter          (simpler alternative for well-separated modes)
        ▼
  TransientSample.q               (scalar modal amplitude signal)
        │
        ▼ find_peaks               (extract local maxima — phase fixing)
        │ discard_initial_transient (remove first n_periods of fast-decaying modes)
        ▼
  Peak sequence (t_peaks, r_peaks)
        │
        ▼ estimate_recovery_rate   (FiniteDifference or EnvelopeODE)
        ▼
  RecoveryRateCurve per sample
        │
        ▼ fit_polynomial_surface   (least-squares polynomial in parameter space)
        ▼
  PolynomialSurface λ(μ, r)
        │
        ▼ extract_diagram          (solve λ = 0 analytically at each amplitude)
        ▼
  BifurcationDiagram
```

The [`forecast`](@ref) function collapses the middle three steps into one call.

## Features

| Feature | API |
|---------|-----|
| ERA-based modal basis extraction | `era_modal_basis` |
| Biorthogonal modal projection | `modal_project` |
| FFT bandpass filter | `bandpass_filter` |
| Peak extraction and phase fixing | `find_peaks`, `discard_initial_transient` |
| Finite-difference recovery-rate estimator | `estimate_recovery_rate(…, FiniteDifference())` |
| Envelope-ODE (Levenberg–Marquardt) estimator | `estimate_recovery_rate(…, EnvelopeODE())` |
| 1st- and 2nd-order polynomial surface fitting | `fit_polynomial_surface` |
| Surface evaluation at new parameter points | `evaluate` |
| Bifurcation diagram extraction | `extract_diagram` |
| Flutter-boundary tracing | `trace_flutter_boundary` |
| Bifurcation type classification | `BifurcationDiagram.bifurcation_type` |
| Model-based state velocity method | `forecast_state_velocity` |
| Latin Hypercube Sampling | `lhs_samples` |
| Amplitude grid generation | `amplitude_grid` |
| Polynomial order convergence check | `convergence_check` |
| Design-matrix condition check | `check_condition` |
| High-level one-call pipeline | `forecast` |

## Running the tests

```julia
using Pkg
Pkg.test("BifurcationForecasting")
```

The test suite covers:

1. **1-DOF analytical** — supercritical and subcritical Hopf normal forms with
   exact solutions; verifies recovery-rate accuracy, surface fitting, flutter
   boundary prediction, and bifurcation-type classification.
2. **4-state aeroelastic typical section** — the full two-parameter test case
   from Riso, Cesnik & Epureanu (2022), testing ERA, modal projection,
   `EnvelopeODE`, and `forecast_state_velocity`.

## Building the documentation

```bash
julia --project=docs docs/make.jl
```

The generated HTML is written to `docs/build/`.

## Licence

Copyright © 2025 David A. W. Barton, University of Bristol.
