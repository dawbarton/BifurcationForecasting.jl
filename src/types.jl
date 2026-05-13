# types.jl
#
# Core data structures for BifurcationForecasting.jl.
#
# Data flow:
#   TransientDataset  →  PolynomialSurface  →  BifurcationDiagram
#
# The user is responsible for any modal projection or channel selection before
# constructing TransientSample; the scalar field q holds the pre-processed
# amplitude signal that enters the recovery-rate pipeline.

# ── Input data ────────────────────────────────────────────────────────────────

"""
    TransientSample(μ, t, q)

One free-decay transient at a fixed control-parameter point.

# Fields
- `μ`: control-parameter vector (length Np).
- `t`: time vector (length N).
- `q`: scalar amplitude signal (length N); must already be modal-projected or
  bandpass-filtered so that it approximates the envelope of the bifurcating mode.
"""
struct TransientSample
    μ::Vector{Float64}   # control-parameter vector (Np,)
    t::Vector{Float64}   # time points (N,)
    q::Vector{Float64}   # scalar signal (N,)
end

"""
    TransientDataset(samples, param_names)

Collection of transient samples at distinct control-parameter values.

# Fields
- `samples`: vector of `TransientSample`.
- `param_names`: names of each control parameter (length Np), used for display.
"""
struct TransientDataset
    samples::Vector{TransientSample}
    param_names::Vector{String}
end

# ── Modal basis (ERA output) ───────────────────────────────────────────────────

"""
    ModalBasis(Φ, Ψ, ΨᵀΦ_inv, eigenvalue)

Bifurcating-mode basis extracted by ERA or direct eigendecomposition.

The real and imaginary parts of the complex right/left eigenvectors are stored
as the two columns of Φ and Ψ respectively.  Pre-computing (ΨᵀΦ)⁻¹ avoids
repeated solves during projection.

# Fields
- `Φ`: (N × 2) matrix `[Re(ϕ_c)  Im(ϕ_c)]`.
- `Ψ`: (N × 2) matrix `[Re(ψ_c)  Im(ψ_c)]`.
- `ΨᵀΦ_inv`: (2 × 2) inverse of `ΨᵀΦ`, pre-computed.
- `eigenvalue`: complex eigenvalue σ_c = g_c + iω_c of the bifurcating mode.

# Reference
Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009, §3.
"""
struct ModalBasis
    Φ::Matrix{Float64}        # (N × 2)
    Ψ::Matrix{Float64}        # (N × 2)
    ΨᵀΦ_inv::Matrix{Float64}  # (2 × 2)
    eigenvalue::ComplexF64    # σ_c = g_c + iω_c
end

# ── Recovery-rate estimation methods ──────────────────────────────────────────

"""
Abstract type for recovery-rate estimation strategies.
Concrete subtypes: `FiniteDifference`, `EnvelopeODE`.
"""
abstract type RecoveryRateMethod end

"""
    FiniteDifference()

Estimate recovery rates from consecutive peak pairs via finite differences of
`log(r)`.  Suitable when many peaks are available (> ~10 per transient).

λ_k = (ln r_{k+1} − ln r_k) / (t_{k+1} − t_k)

# Reference
Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009, §4.2 (Eq. 9).
"""
struct FiniteDifference <: RecoveryRateMethod end

"""
    EnvelopeODE(poly_order)

Estimate recovery rates by fitting the amplitude-envelope ODE

    ṙ = r (λ₀ + λ₁r + … + λₚrᵖ)

to observed peak amplitudes via Levenberg–Marquardt nonlinear least squares
(NonlinearSolve.jl).  Preferred for oscillatory systems with few peaks.

# Fields
- `poly_order`: polynomial degree p (default 4; increase until results stabilise).

# Reference
Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009, §4.3 (Eq. 11).
"""
struct EnvelopeODE <: RecoveryRateMethod
    poly_order::Int
end
EnvelopeODE() = EnvelopeODE(4)

# ── Recovery-rate data ────────────────────────────────────────────────────────

"""
    RecoveryRateCurve(μ, r, λ, method)

Recovery rate λ(r) estimated at a single control-parameter point μ.

After estimation (either finite-difference or envelope-ODE fitting), the curve
is represented as paired (r, λ) vectors that are subsequently interpolated onto
a common amplitude grid by `fit_polynomial_surface`.

# Fields
- `μ`: control-parameter vector for this sample.
- `r`: amplitude values (not necessarily uniform).
- `λ`: recovery rates at those amplitudes.
- `method`: the `RecoveryRateMethod` used to produce this curve.
"""
struct RecoveryRateCurve
    μ::Vector{Float64}
    r::Vector{Float64}
    λ::Vector{Float64}
    method::RecoveryRateMethod
end

# ── Polynomial surface (primary forecast output) ──────────────────────────────

"""
    PolynomialSurface(r_grid, coeffs, param_mean, param_scale,
                      order, Np, Nc, flutter_frequency)

Fitted polynomial approximation of the recovery-rate surface λ(μ, r).

At each amplitude r̃_k on the grid, the recovery rate is approximated as a
polynomial in the (normalised) control parameters:

    λ(r̃_k, μ) ≈ cᵀ φ(μ̂)

where μ̂ = (μ − param_mean) / param_scale and φ is the monomial basis vector.

**Coefficient ordering** (matches `_build_design_row`):
- order = 1: `[a₀, a₁⁽¹⁾, …, a₁⁽ᴺᵖ⁾]`          (Nc = 1 + Np)
- order = 2: `[a₀, a₁⁽¹⁾, …, a₁⁽ᴺᵖ⁾, a₂⁽¹¹⁾, a₂⁽¹²⁾, …, a₂⁽ᴺᵖᴺᵖ⁾]`
             where the a₂ entries follow the upper-triangle of the symmetric
             matrix in row-major order  (Nc = 1 + Np + Np(Np+1)/2).

# Fields
- `r_grid`: amplitude grid (Nr,).
- `coeffs`: (Nc × Nr) matrix of polynomial coefficients per amplitude slice.
- `param_mean`: mean of parameter samples used for normalisation (Np,).
- `param_scale`: std of parameter samples used for normalisation (Np,).
- `order`: Taylor expansion order in parameter space (1 or 2).
- `Np`: number of control parameters.
- `Nc`: number of polynomial coefficients per amplitude slice.
- `flutter_frequency`: ω_c (rad/time), by-product of ERA or eigenvalue step.

# Reference
Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §2.
"""
struct PolynomialSurface
    r_grid::Vector{Float64}       # (Nr,)
    coeffs::Matrix{Float64}       # (Nc × Nr)
    param_mean::Vector{Float64}   # (Np,)
    param_scale::Vector{Float64}  # (Np,)
    order::Int
    Np::Int
    Nc::Int
    flutter_frequency::Float64
end

# ── Bifurcation diagram ───────────────────────────────────────────────────────

"""
    BifurcationPoint(μ, r, λ_slope)

A single point on the forecasted bifurcation diagram.

# Fields
- `μ`: control-parameter vector at this point.
- `r`: amplitude (0 for the flutter boundary).
- `λ_slope`: ∂λ/∂r evaluated at (μ, r) along the zero-level set.
  Negative ⟹ stable LCO; positive ⟹ unstable LCO.
"""
struct BifurcationPoint
    μ::Vector{Float64}
    r::Float64
    λ_slope::Float64
end

"""
    BifurcationDiagram(points, flutter_boundary, bifurcation_type, flutter_frequency)

Forecasted bifurcation diagram in (μ, r) space.

# Fields
- `points`: vector of `BifurcationPoint` covering all amplitude levels.
- `flutter_boundary`: control-parameter vectors where r = 0 (flutter onset).
- `bifurcation_type`: `:supercritical`, `:subcritical`, or `:unknown`.
- `flutter_frequency`: ω_c (rad/time) of the bifurcating mode.

# Reference
Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §2.3.
"""
struct BifurcationDiagram
    points::Vector{BifurcationPoint}
    flutter_boundary::Vector{Vector{Float64}}
    bifurcation_type::Symbol
    flutter_frequency::Float64
end
