# diagram.jl
#
# Bifurcation diagram extraction from a fitted PolynomialSurface.
#
# For each amplitude r̃_k on the grid, the zero-recovery-rate condition
# λ(μ, r̃_k) = 0 is solved in the normalised parameter space.  With a
# scalar parameter this is a quadratic in μ̂; with multiple parameters the
# search is parameterised as μ̂ = μ̂★ + δ̃ μ̄ and reduces to a scalar
# quadratic in δ̃ (Riso, Cesnik & Epureanu 2021, §2.3).
#
# Public API
# ──────────
#   extract_diagram(surface; μ_ref, μ_dir, r_min) → BifurcationDiagram
#   trace_flutter_boundary(surface; μ_ref, μ_dir) → Vector{Vector{Float64}}
#
# Private helpers
# ───────────────
#   _solve_quadratic_δ  — solves the scalar quadratic in δ̃
#   _bifurcation_type   — classifies super/subcritical from LCO position relative to flutter boundary
#
# References
# ──────────
#   Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §2.3
#     https://doi.org/10.1016/j.jfluidstructs.2020.103201
#   bifurcation_forecasting.md §4.6–4.7, §9

# ── Public: full diagram extraction ──────────────────────────────────────────

"""
    extract_diagram(surface; μ_ref, μ_dir, r_min) → BifurcationDiagram

Extract the bifurcation diagram from a fitted `PolynomialSurface` by finding
the zero-recovery-rate locus λ(μ, r) = 0 at each amplitude grid point.

The zero-crossing condition is solved parametrically.  A reference parameter
point `μ_ref` and unit direction `μ_dir` define the search line

    μ̂(δ̃) = μ̂★ + δ̃ μ̄

in normalised parameter space.  The scalar δ̃ is found by solving the
resulting quadratic analytically.

# Keyword arguments
- `μ_ref`: reference parameter vector (unnormalised) for the search origin.
  Default: the mean of the parameter samples used for fitting.
- `μ_dir`: search direction in parameter space (unnormalised unit vector).
  Default: first standard basis vector e₁.
- `r_min`: minimum amplitude to include (default 0.0).

# Returns
A `BifurcationDiagram` containing:
- `points`: one `BifurcationPoint` per grid amplitude where a real root exists.
- `flutter_boundary`: parameter values where r → 0 (the flutter boundary).
- `bifurcation_type`: `:supercritical`, `:subcritical`, or `:unknown`.
- `flutter_frequency`: ω_c stored in the surface.
"""
function extract_diagram(
        surface::PolynomialSurface;
        μ_ref::Union{Nothing, AbstractVector{<:Real}} = nothing,
        μ_dir::Union{Nothing, AbstractVector{<:Real}} = nothing,
        r_min::Real = 0.0
    )::BifurcationDiagram

    Np = surface.Np
    Nr = length(surface.r_grid)

    # Default reference point: parameter mean (zero in normalised space).
    μ★_norm = if isnothing(μ_ref)
        zeros(Np)
    else
        (collect(Float64, μ_ref) .- surface.param_mean) ./ surface.param_scale
    end

    # Default search direction: e₁ in normalised space.
    μ̄_norm = if isnothing(μ_dir)
        e1 = zeros(Np); e1[1] = 1.0; e1
    else
        d = collect(Float64, μ_dir) ./ surface.param_scale
        d ./ norm(d)
    end

    points = BifurcationPoint[]

    for k in 1:Nr
        r_k = surface.r_grid[k]
        r_k < r_min && continue

        c = surface.coeffs[:, k]   # coefficient vector (Nc,)
        a0 = c[1]
        a1 = c[2:(1 + Np)]
        A2 = extract_a2(surface, k)   # (Np × Np) in normalised space

        # Evaluate polynomial at μ★ to get the constant term of the quadratic.
        c_const = a0 + dot(a1, μ★_norm) + dot(μ★_norm, A2 * μ★_norm)
        c_lin = dot(a1, μ̄_norm) + 2.0 * dot(μ★_norm, A2 * μ̄_norm)
        c_quad = dot(μ̄_norm, A2 * μ̄_norm)

        δ̃ = _solve_quadratic_δ(c_quad, c_lin, c_const, surface.order)
        isnan(δ̃) && continue

        # Reconstruct the bifurcation point in unnormalised parameter space.
        μ̂_bif = μ★_norm .+ δ̃ .* μ̄_norm
        μ_bif = μ̂_bif .* surface.param_scale .+ surface.param_mean

        # Compute ∂λ/∂r at the bifurcation point to assess LCO stability.
        # Use a two-point finite difference across amplitude indices where possible.
        λ_slope = if k < Nr
            r_next = surface.r_grid[k + 1]
            c_next = surface.coeffs[:, k + 1]
            a0n = c_next[1]; a1n = c_next[2:(1 + Np)]; A2n = extract_a2(surface, k + 1)
            λ_next = a0n + dot(a1n, μ̂_bif) + dot(μ̂_bif, A2n * μ̂_bif)
            λ_next / (r_next - r_k)   # ∂λ/∂r ≈ (λ(r+Δr) - 0) / Δr
        else
            0.0
        end

        push!(points, BifurcationPoint(μ_bif, r_k, λ_slope))
    end

    flutter_boundary = trace_flutter_boundary(surface; μ_ref = μ_ref, μ_dir = μ_dir)

    # Search direction (unnormalised) for bifurcation type comparison.
    μ_dir_eff = if isnothing(μ_dir)
        e1 = zeros(Np); e1[1] = 1.0; e1
    else
        collect(Float64, μ_dir)
    end
    bif_type = _bifurcation_type(points, surface, flutter_boundary, μ_dir_eff)

    return BifurcationDiagram(
        points, flutter_boundary, bif_type,
        surface.flutter_frequency
    )
end

# ── Public: flutter boundary tracing ─────────────────────────────────────────

"""
    trace_flutter_boundary(surface; μ_ref, μ_dir) → Vector{Vector{Float64}}

Find parameter values on the flutter boundary (r = 0 zero of λ).

For a single-parameter surface the boundary is a point; for multi-parameter
surfaces this returns the zero contour along the specified search direction.
Uses the r = 0 limit, i.e. only the a₀ and a₁ terms (the quadratic a₂ term
vanishes at r = 0 in the small-amplitude limit and is numerically noisy).

Returns a vector of unnormalised parameter vectors on the flutter boundary.
"""
function trace_flutter_boundary(
        surface::PolynomialSurface;
        μ_ref::Union{Nothing, AbstractVector{<:Real}} = nothing,
        μ_dir::Union{Nothing, AbstractVector{<:Real}} = nothing
    )::Vector{Vector{Float64}}

    Np = surface.Np

    μ★_norm = if isnothing(μ_ref)
        zeros(Np)
    else
        (collect(Float64, μ_ref) .- surface.param_mean) ./ surface.param_scale
    end

    μ̄_norm = if isnothing(μ_dir)
        e1 = zeros(Np); e1[1] = 1.0; e1
    else
        d = collect(Float64, μ_dir) ./ surface.param_scale
        d ./ norm(d)
    end

    # Evaluate at r = 0 using only the first amplitude grid slice (nearest to zero).
    # The flutter boundary condition is a₀(0) + a₁(0)ᵀμ = 0.  The quadratic term
    # is included if order == 2 for consistency, but this is less reliable for
    # the boundary (see §8 of bifurcation_forecasting.md).
    k = 1   # use the smallest amplitude on the grid
    c = surface.coeffs[:, k]
    a0 = c[1]
    a1 = c[2:(1 + Np)]
    A2 = extract_a2(surface, k)

    c_const = a0 + dot(a1, μ★_norm) + dot(μ★_norm, A2 * μ★_norm)
    c_lin = dot(a1, μ̄_norm) + 2.0 * dot(μ★_norm, A2 * μ̄_norm)
    c_quad = dot(μ̄_norm, A2 * μ̄_norm)

    δ̃ = _solve_quadratic_δ(c_quad, c_lin, c_const, surface.order)
    isnan(δ̃) && return Vector{Float64}[]

    μ̂_fb = μ★_norm .+ δ̃ .* μ̄_norm
    μ_fb = μ̂_fb .* surface.param_scale .+ surface.param_mean

    return [μ_fb]
end

# ── Private helpers ───────────────────────────────────────────────────────────

# Solve the scalar quadratic (or linear) equation in δ̃:
#   c_quad * δ̃² + c_lin * δ̃ + c_const = 0
# When order == 1, c_quad is ignored and the equation is linear.
# Returns NaN if no real root exists or the linear equation has no solution.
# When two real roots exist, the one with smaller |δ̃| is returned (closest to μ★).
function _solve_quadratic_δ(
        c_quad::Real, c_lin::Real, c_const::Real,
        order::Int
    )::Float64

    if order == 1 || c_quad == 0.0
        # Linear equation: c_lin * δ̃ + c_const = 0
        abs(c_lin) < 1.0e-14 && return NaN
        return -c_const / c_lin
    else
        # Quadratic formula.
        disc = c_lin^2 - 4.0 * c_quad * c_const
        disc < 0.0 && return NaN
        sq = sqrt(disc)
        δ1 = (-c_lin + sq) / (2.0 * c_quad)
        δ2 = (-c_lin - sq) / (2.0 * c_quad)
        # Return the root closest to the reference (smallest |δ̃|), i.e. the
        # physically nearest bifurcation point to the measured regime.
        return abs(δ1) <= abs(δ2) ? δ1 : δ2
    end
end

# Classify the bifurcation type from the slope of λ with respect to r at the
# diagram points.  Supercritical: λ slope negative (λ decreases as r grows,
# meaning the LCO is born stable).  Subcritical: λ slope positive near r = 0
# (unstable LCO branch exists below flutter onset).
function _bifurcation_type(
        points::Vector{BifurcationPoint},
        surface::PolynomialSurface,
        flutter_boundary::Vector{Vector{Float64}},
        μ_dir::AbstractVector{<:Real}
    )::Symbol
    isempty(points) && return :unknown

    # Primary criterion: majority vote comparing LCO branch position to flutter boundary
    # along μ_dir over the whole r_grid.  Supercritical LCOs predominantly sit above
    # the flutter point; subcritical ones below.  Using all points avoids the noise
    # that dominates ∂λ/∂r at very small r (where the nonlinear correction ~r² is
    # smaller than the surface fit residual).
    if !isempty(flutter_boundary)
        μ_fb = flutter_boundary[1]
        dir = μ_dir ./ norm(μ_dir)
        δ_bif = [dot(p.μ .- μ_fb, dir) for p in points]
        n_pos = count(>(0.0), δ_bif)
        n_neg = count(<(0.0), δ_bif)
        n_pos > n_neg && return :supercritical
        n_neg > n_pos && return :subcritical
    end

    # Fallback: majority vote on λ_slope sign across all points.
    n_neg_slope = count(p -> p.λ_slope < 0.0, points)
    n_pos_slope = count(p -> p.λ_slope > 0.0, points)
    n_neg_slope > n_pos_slope && return :supercritical
    n_pos_slope > n_neg_slope && return :subcritical
    return :unknown
end
