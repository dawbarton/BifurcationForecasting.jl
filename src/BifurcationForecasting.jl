# BifurcationForecasting.jl
#
# Bifurcation forecasting for flutter and limit-cycle oscillation (LCO) onset,
# based on the critical-slowing-down (CSD) method of the Epureanu group.
#
# The core idea: as a control parameter approaches the flutter boundary, a
# system's free-decay transients slow down (CSD).  Quantifying this slow-down
# across a small set of pre-flutter parameter samples allows the bifurcation
# diagram to be reconstructed by extrapolation, without ever crossing into the
# unstable regime and without requiring a system model (transient-data method).
#
# A model-based variant (state velocity method) is also provided.
#
# Primary references
# ──────────────────
#   Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009
#     https://doi.org/10.1115/1.4033920
#   Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201
#     https://doi.org/10.1016/j.jfluidstructs.2020.103201  [multi-parameter]
#   Riso, Cesnik & Epureanu (2022) AIAA J. 60, 5401–5413
#     https://doi.org/10.2514/1.J061860  [state velocity method]

module BifurcationForecasting

using LinearAlgebra
using Random
using Statistics
using FFTW
using Interpolations
using OrdinaryDiffEq
using NonlinearSolve

include("types.jl")
include("era.jl")
include("signal.jl")
include("recovery_rate.jl")
include("surface.jl")
include("diagram.jl")
include("state_velocity.jl")
include("sampling.jl")
include("diagnostics.jl")
include("pipeline.jl")

# ── Public API ────────────────────────────────────────────────────────────────

# Types
export TransientSample, TransientDataset
export ModalBasis
export FiniteDifference, EnvelopeODE
export RecoveryRateCurve
export PolynomialSurface
export BifurcationPoint, BifurcationDiagram

# ERA and modal filtering
export era_modal_basis, modal_project, bandpass_filter

# Signal processing
export find_peaks, discard_initial_transient

# Recovery rate estimation
export estimate_recovery_rate

# Surface fitting
export fit_polynomial_surface, evaluate, extract_a2

# Bifurcation diagram extraction
export extract_diagram, trace_flutter_boundary

# State velocity method
export forecast_state_velocity

# Sampling utilities
export lhs_samples, amplitude_grid

# Diagnostics
export convergence_check, check_condition

# High-level pipeline
export forecast

end # module BifurcationForecasting
