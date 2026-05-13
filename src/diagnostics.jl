# diagnostics.jl
#
# Diagnostic utilities for the bifurcation forecasting pipeline.
#
# Two checks are provided:
#
#   convergence_check — compare first- vs. second-order surface fits to decide
#     whether a quadratic Taylor expansion is necessary.
#
#   check_condition   — report the condition number of the design matrix for
#     each amplitude slice of a PolynomialSurface, flagging ill-conditioning.
#
# References
# ──────────
#   Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §8
#     https://doi.org/10.1016/j.jfluidstructs.2020.103201  (convergence check)
#   bifurcation_forecasting.md §8 (Polynomial Fitting Stability)

"""
    convergence_check(curves, r_grid; tol, μ_ref, μ_dir) → NamedTuple

Compare first- and second-order polynomial surface fits to assess whether the
quadratic expansion is necessary.

The check (Riso 2021 §8 / algorithm doc §8) computes the relative difference
in the predicted flutter speed (the root of λ = 0 at r = 0) between the two
fits.  If the difference is below `tol`, the first-order fit is deemed
sufficient.

# Arguments
- `curves`: vector of `RecoveryRateCurve` (as returned by the recovery-rate
  estimation step).
- `r_grid`: amplitude grid for the surface fit.

# Keyword arguments
- `tol`: relative tolerance for declaring convergence (default 0.01, i.e. 1%).
- `μ_ref`, `μ_dir`: reference point and search direction for diagram extraction
  (passed to `extract_diagram`; see that function's documentation).

# Returns
A `NamedTuple` with fields:
- `converged::Bool` — whether first-order is sufficient.
- `rel_diff::Float64` — relative difference in predicted flutter speed.
- `surface_1::PolynomialSurface` — first-order surface.
- `surface_2::PolynomialSurface` — second-order surface.
- `diagram_1::BifurcationDiagram` — bifurcation diagram from first-order fit.
- `diagram_2::BifurcationDiagram` — bifurcation diagram from second-order fit.
"""
function convergence_check(
        curves::AbstractVector{RecoveryRateCurve},
        r_grid::AbstractVector{<:Real};
        tol::Real = 0.01,
        μ_ref::Union{Nothing, AbstractVector{<:Real}} = nothing,
        μ_dir::Union{Nothing, AbstractVector{<:Real}} = nothing
    )

    surf1 = fit_polynomial_surface(curves, r_grid; order = 1)
    surf2 = fit_polynomial_surface(curves, r_grid; order = 2)

    diag1 = extract_diagram(surf1; μ_ref, μ_dir)
    diag2 = extract_diagram(surf2; μ_ref, μ_dir)

    # Extract flutter boundary parameters from each diagram; use the first point.
    μ_fb1 = isempty(diag1.flutter_boundary) ? NaN : diag1.flutter_boundary[1][1]
    μ_fb2 = isempty(diag2.flutter_boundary) ? NaN : diag2.flutter_boundary[1][1]

    rel_diff = abs(μ_fb1 - μ_fb2) / max(abs(μ_fb1), 1.0e-12)
    converged = rel_diff < tol

    return (;
        converged, rel_diff, surface_1 = surf1, surface_2 = surf2,
        diagram_1 = diag1, diagram_2 = diag2,
    )
end

"""
    check_condition(curves, r_grid; order, warn_threshold) → Vector{Float64}

Compute the condition number of the polynomial design matrix for each amplitude
slice and warn if any exceed `warn_threshold`.

A large condition number (> ~1e6) indicates that the parameter samples are
poorly conditioned for the chosen polynomial order; consider adding more
samples, reducing the order, or enabling Tikhonov regularisation.

# Arguments
- `curves`, `r_grid`, `order`: same meaning as in `fit_polynomial_surface`.

# Keyword arguments
- `warn_threshold`: condition number above which a warning is printed
  (default `1e6`).

# Returns
Vector of condition numbers (length Nr, one per amplitude slice).
"""
function check_condition(
        curves::AbstractVector{RecoveryRateCurve},
        r_grid::AbstractVector{<:Real};
        order::Int = 1,
        warn_threshold::Real = 1.0e6
    )::Vector{Float64}

    Ns = length(curves)
    Np = length(curves[1].μ)

    # Build the design matrix (same normalisation as fit_polynomial_surface).
    μ_mat = stack([c.μ for c in curves])
    param_mean = vec(mean(μ_mat, dims = 2))
    param_scale = vec(std(μ_mat, dims = 2))
    param_scale[param_scale .== 0.0] .= 1.0

    A = zeros(Ns, order == 1 ? 1 + Np : 1 + Np + Np * (Np + 1) ÷ 2)
    for l in 1:Ns
        μ_norm = (curves[l].μ .- param_mean) ./ param_scale
        A[l, :] = _build_design_row(μ_norm, order)
    end

    κ = cond(A)   # condition number of A (same for all amplitude slices)
    Nr = length(r_grid)

    if κ > warn_threshold
        @warn "check_condition: design matrix condition number $κ exceeds threshold $warn_threshold. " *
            "Consider adding samples, reducing polynomial order, or enabling ridge regularisation."
    end

    return fill(κ, Nr)
end
