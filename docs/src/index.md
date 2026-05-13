# BifurcationForecasting.jl

!!! warning "AI-generated code"
    This package was written by [Claude](https://www.anthropic.com) (Anthropic)
    under the direction of David A. W. Barton.  All scientific content has been
    reviewed by the author, but the implementation should be treated with
    appropriate scepticism and validated independently before use in research.

**BifurcationForecasting.jl** implements the *critical-slowing-down* (CSD)
method for predicting flutter onset and limit-cycle oscillation (LCO) amplitude
without crossing the stability boundary.

## Background

As a control parameter (e.g. airspeed) approaches a Hopf bifurcation, the
dominant eigenvalue of the linearised system approaches the imaginary axis.
Free-decay transients initiated in the pre-flutter regime therefore slow down —
a phenomenon known as *critical slowing down*.  The recovery rate

```
λ(μ, r) = ṙ / r
```

quantifies how quickly the system returns to equilibrium from amplitude `r`
at parameter value `μ`.  On the bifurcation diagram, `λ = 0` defines:

- **Flutter boundary** — where `r = 0` (linear instability onset).
- **LCO branch** — where `r > 0` (finite-amplitude steady oscillation).

Measuring `λ(r)` curves from a small number of pre-flutter transients allows
the zero-level set to be reconstructed by extrapolation, giving a
*data-driven bifurcation forecast* without a system model and without
destabilising the system.

A *model-based* variant ([`forecast_state_velocity`](@ref)) is also provided
for cases where the nonlinear state equations are known but time integration
is expensive.

## Primary references

| Reference | Topic |
|:----------|:------|
| Ghadami & Epureanu (2016) *J. Comput. Nonlinear Dyn.* **11**, 061009. [doi](https://doi.org/10.1115/1.4033920) | Core CSD method, ERA projection, finite-difference and ODE envelope estimators |
| Riso, Cesnik & Epureanu (2021) *J. Fluids Struct.* **101**, 103201. [doi](https://doi.org/10.1016/j.jfluidstructs.2020.103201) | Multi-parameter polynomial surface, bifurcation diagram extraction |
| Riso, Cesnik & Epureanu (2022) *AIAA J.* **60**, 5401–5413. [doi](https://doi.org/10.2514/1.J061860) | State velocity method |
| Juang & Pappa (1985) *J. Guid. Control Dyn.* **8**, 620–627 | Eigensystem Realisation Algorithm (ERA) |

## Package contents

```
BifurcationForecasting.jl
├── types.jl          — core data structures
├── era.jl            — ERA, modal projection, bandpass filtering
├── signal.jl         — peak extraction, initial-transient discard
├── recovery_rate.jl  — FiniteDifference and EnvelopeODE estimators
├── surface.jl        — polynomial surface fitting and evaluation
├── diagram.jl        — bifurcation diagram extraction
├── state_velocity.jl — model-based (state velocity) pipeline
├── sampling.jl       — LHS and amplitude-grid utilities
├── diagnostics.jl    — convergence check, condition-number check
└── pipeline.jl       — high-level forecast() entry point
```

## Installation

The package is not registered.  From the Julia REPL:

```julia
using Pkg
Pkg.add(url = "https://github.com/dawbarton/BifurcationForecasting.jl")
```

## Quick start

```julia
using BifurcationForecasting

# 1. Wrap pre-processed transients (one per parameter point).
samples = [TransientSample([μ], t_vec, q_vec) for (μ, t_vec, q_vec) in data]
dataset = TransientDataset(samples, ["airspeed"])

# 2. Run the full pipeline (requires flutter-mode frequency ω_c).
surface = forecast(dataset, ω_c)

# 3. Extract the bifurcation diagram.
diagram = extract_diagram(surface)
println("Flutter boundary: ", diagram.flutter_boundary)
println("Bifurcation type: ", diagram.bifurcation_type)
```

See the [Workflow](@ref) page for a self-contained worked example.
