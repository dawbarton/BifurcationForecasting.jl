# surface.jl
#
# Fitting the recovery-rate polynomial surface λ(μ, r) from a collection of
# RecoveryRateCurve objects (one per control-parameter sample).
#
# The surface is the core intermediate output of the bifurcation forecasting
# pipeline.  At each amplitude grid point r̃_k the recovery rate is
# approximated as a polynomial in the (normalised) control-parameter vector:
#
#   λ(r̃_k, μ) ≈ a₀ + a₁ᵀμ̂ + μ̂ᵀ a₂ μ̂
#
# where μ̂ = (μ − mean(μ)) / std(μ) and a₂ is symmetric.  All coefficients
# are fitted by ordinary least squares (via Julia's built-in QR path `\`).
#
# Public API
# ──────────
#   fit_polynomial_surface(curves, r_grid; order, ridge) → PolynomialSurface
#   evaluate(surface, μ) → Vector{Float64}  (λ at every r in r_grid)
#   extract_a2(surface, k) → Matrix{Float64}  (symmetric a₂ at grid index k)
#
# Private helpers
# ───────────────
#   _build_design_row  — monomial row of the design matrix for one sample
#   _interp_curve      — interpolate a RecoveryRateCurve onto an r grid
#   _assemble_λ_matrix — collect interpolated λ values into (Ns × Nr) matrix
#
# References
# ──────────
#   Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §2
#     https://doi.org/10.1016/j.jfluidstructs.2020.103201
#   bifurcation_forecasting.md §4.4–4.5

# ── Public: surface fitting ───────────────────────────────────────────────────

"""
    fit_polynomial_surface(curves, r_grid; order, ridge) → PolynomialSurface

Fit the recovery-rate polynomial surface λ(μ, r) from a collection of
per-sample recovery-rate curves.

At each amplitude r̃_k in `r_grid`, the λ values from all samples are
interpolated and then fitted by the polynomial

    λ(r̃_k, μ) ≈ a₀ + a₁ᵀμ̂ + μ̂ᵀ a₂ μ̂          (order 2)
    λ(r̃_k, μ) ≈ a₀ + a₁ᵀμ̂                         (order 1)

where μ̂ = (μ − mean(μ)) / std(μ) is the normalised parameter vector.
Fitting is done by least squares; `ridge > 0` adds Tikhonov regularisation
(`ridge × I`) to stabilise ill-conditioned design matrices.

# Arguments
- `curves`: vector of `RecoveryRateCurve` (length Ns, one per parameter sample).
- `r_grid`: amplitude values at which to evaluate and store the surface (Nr,).

# Keyword arguments
- `order`: Taylor expansion order in parameter space (1 or 2; default 1).
- `ridge`: Tikhonov regularisation parameter (default 0.0).

# Returns
A `PolynomialSurface` with `coeffs` of size `(Nc × Nr)`.

# Reference
Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §2.
"""
function fit_polynomial_surface(
        curves::AbstractVector{RecoveryRateCurve},
        r_grid::AbstractVector{<:Real};
        order::Int = 1,
        ridge::Real = 0.0
    )::PolynomialSurface

    Ns = length(curves)
    Nr = length(r_grid)
    Ns < 1 && error("fit_polynomial_surface: need at least one curve; got $Ns")
    (order == 1 || order == 2) || error("order must be 1 or 2; got $order")

    # All curves must have the same parameter dimension.
    Np = length(curves[1].μ)
    all(length(c.μ) == Np for c in curves) ||
        error("all RecoveryRateCurve.μ must have the same length")

    Nc = order == 1 ? 1 + Np : 1 + Np + Np * (Np + 1) ÷ 2

    # Normalise parameters: μ̂ = (μ − mean) / std across all samples.
    μ_mat = stack([c.μ for c in curves])   # (Np × Ns)
    param_mean = vec(mean(μ_mat, dims = 2))
    param_scale = vec(std(μ_mat, dims = 2))
    # Avoid division by zero when a parameter is constant across samples.
    param_scale[param_scale .== 0.0] .= 1.0

    # Assemble design matrix A (Ns × Nc): one design row per sample.
    A = zeros(Ns, Nc)
    for l in 1:Ns
        μ_norm = (curves[l].μ .- param_mean) ./ param_scale
        A[l, :] = _build_design_row(μ_norm, order)
    end

    # Interpolate all curves onto r_grid: λ_mat is (Ns × Nr).
    λ_mat = _assemble_λ_matrix(curves, collect(Float64, r_grid))

    # Fit coefficients for each amplitude slice by least squares (with optional ridge).
    # Stores result as (Nc × Nr); column k holds coefficients for r_grid[k].
    coeffs = zeros(Nc, Nr)
    if ridge > 0.0
        AᵀA = A' * A + ridge * I(Nc)
        Aᵀ = A'
        for k in 1:Nr
            coeffs[:, k] = AᵀA \ (Aᵀ * λ_mat[:, k])
        end
    else
        # pinv handles rank-deficient cases (e.g. constant parameter across all samples)
        # without error; gives the minimum-norm least-squares solution.
        coeffs .= pinv(A) * λ_mat
    end

    r_grid_f = collect(Float64, r_grid)

    # Extract flutter frequency from the first curve's modal basis eigenvalue if
    # available; otherwise NaN.  The pipeline sets this on PolynomialSurface later.
    flutter_frequency = NaN

    return PolynomialSurface(
        r_grid_f, coeffs, param_mean, param_scale,
        order, Np, Nc, flutter_frequency
    )
end

# ── Public: evaluation ────────────────────────────────────────────────────────

"""
    evaluate(surface::PolynomialSurface, μ) → Vector{Float64}

Evaluate the fitted recovery-rate surface at a new control-parameter point μ.

Returns a vector of recovery rates at each amplitude in `surface.r_grid`.

# Example
```julia
λ_vals = evaluate(surf, [0.8])  # λ(r) at each r in surf.r_grid for μ = 0.8
```
"""
function evaluate(
        surface::PolynomialSurface,
        μ::AbstractVector{<:Real}
    )::Vector{Float64}
    length(μ) == surface.Np ||
        error("μ has length $(length(μ)); surface expects Np=$(surface.Np)")
    μ_norm = (μ .- surface.param_mean) ./ surface.param_scale
    row = _build_design_row(μ_norm, surface.order)
    # coeffs is (Nc × Nr); return (Nr,) vector.
    return vec(row' * surface.coeffs)
end

# ── Public: coefficient extraction ───────────────────────────────────────────

"""
    extract_a2(surface::PolynomialSurface, k::Int) → Matrix{Float64}

Extract the symmetric second-order coefficient matrix a₂ at amplitude grid
index k (1-based).  Returns an (Np × Np) zero matrix if `surface.order == 1`.

The coefficients are in normalised parameter space (μ̂ = (μ − mean)/scale).
Downstream zero-crossing solvers must work in the same normalised space.
"""
function extract_a2(surface::PolynomialSurface, k::Int)::Matrix{Float64}
    Np = surface.Np
    surface.order == 1 && return zeros(Np, Np)
    c = surface.coeffs[:, k]
    A2 = zeros(Np, Np)
    idx = 2 + Np   # first a₂ coefficient starts after [a₀, a₁⁽¹⁾,…,a₁⁽ᴺᵖ⁾]
    for i in 1:Np
        for j in i:Np
            A2[i, j] = c[idx]
            A2[j, i] = c[idx]   # symmetric
            idx += 1
        end
    end
    return A2
end

# ── Private helpers ───────────────────────────────────────────────────────────

# Build one row of the polynomial design matrix for a normalised parameter
# vector μ̂.  Coefficient ordering (matches PolynomialSurface docstring):
#   order=1: [1, μ̂₁, …, μ̂_Np]
#   order=2: [1, μ̂₁, …, μ̂_Np, μ̂₁², μ̂₁μ̂₂, …, μ̂₁μ̂_Np, μ̂₂², …, μ̂_Np²]
#            (upper-triangle in row-major order for the quadratic terms)
function _build_design_row(
        μ_norm::AbstractVector{<:Real},
        order::Int
    )::Vector{Float64}
    Np = length(μ_norm)
    Nc = order == 1 ? 1 + Np : 1 + Np + Np * (Np + 1) ÷ 2
    row = Vector{Float64}(undef, Nc)
    row[1] = 1.0
    row[2:(1 + Np)] = μ_norm
    if order == 2
        idx = 2 + Np
        for i in 1:Np
            for j in i:Np
                row[idx] = μ_norm[i] * μ_norm[j]
                idx += 1
            end
        end
    end
    return row
end

# Interpolate a single RecoveryRateCurve onto the supplied amplitude grid.
# The curve's (r, λ) pairs may not be in ascending order; they are sorted
# before interpolation.  Linear interpolation is used; outside the curve's
# range, the boundary slopes are extrapolated linearly.
function _interp_curve(
        curve::RecoveryRateCurve,
        r_grid::AbstractVector{<:Real}
    )::Vector{Float64}

    r_src = curve.r
    λ_src = curve.λ

    # Sort by ascending r (FiniteDifference curves have decreasing r).
    perm = sortperm(r_src)
    r_srt = r_src[perm]
    λ_srt = λ_src[perm]

    # Remove duplicate r values (keep last; duplicates arise rarely from finite diff).
    keep = [true; diff(r_srt) .> 0]
    r_srt = r_srt[keep]
    λ_srt = λ_srt[keep]

    length(r_srt) < 2 && return fill(NaN, length(r_grid))

    itp = linear_interpolation(r_srt, λ_srt; extrapolation_bc = Line())
    return [itp(r) for r in r_grid]
end

# Assemble the (Ns × Nr) matrix of λ values.
# Row l corresponds to sample l; column k corresponds to r_grid[k].
function _assemble_λ_matrix(
        curves::AbstractVector{RecoveryRateCurve},
        r_grid::AbstractVector{Float64}
    )::Matrix{Float64}
    Ns = length(curves)
    Nr = length(r_grid)
    λ_mat = zeros(Ns, Nr)
    for l in 1:Ns
        λ_mat[l, :] = _interp_curve(curves[l], r_grid)
    end
    return λ_mat
end
