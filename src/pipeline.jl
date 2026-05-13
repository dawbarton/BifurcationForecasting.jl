# pipeline.jl
#
# High-level entry point for the bifurcation forecasting pipeline.
#
# `forecast` wraps the full transient-data method in a single call:
#
#   TransientDataset
#     → modal projection or bandpass filtering (optional, caller's choice)
#     → peak extraction + transient discard (find_peaks, discard_initial_transient)
#     → recovery-rate estimation (FiniteDifference or EnvelopeODE)
#     → polynomial surface fitting (fit_polynomial_surface)
#     → PolynomialSurface
#
# The caller is responsible for pre-processing (ERA or bandpass filtering)
# before constructing TransientSample objects.  This keeps the pipeline
# composable: each step can be invoked individually for diagnostics.
#
# References
# ──────────
#   Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009
#   Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201
#   bifurcation_forecasting.md §5 (complete algorithm pseudocode)

"""
    forecast(dataset::TransientDataset, ω_c;
             method, r_grid, order, ridge,
             min_prominence, n_periods,
             λ_init, r_eval, abstol, reltol) → PolynomialSurface

Run the full transient-data bifurcation forecasting pipeline.

Applies phase fixing, recovery-rate estimation, and polynomial surface fitting
to all samples in `dataset`.  The caller must have already projected each
`TransientSample.q` onto the bifurcating mode (via ERA or bandpass filtering)
before passing to this function.

# Arguments
- `dataset`: collection of pre-processed transient samples.
- `ω_c`: angular frequency of the bifurcating mode (rad/time), used to
  determine the initial-transient discard window.

# Keyword arguments
- `method`: recovery-rate estimation method (default `EnvelopeODE()`).
- `r_grid`: amplitude grid for the surface; if `nothing`, a 50-point linear
  grid from the minimum to maximum peak amplitude across all samples is used.
- `order`: polynomial order for the surface fit (1 or 2; default 1).
- `ridge`: Tikhonov regularisation for the surface fit (default 0.0).
- `min_prominence`: minimum peak height (default 0.0).
- `n_periods`: number of initial oscillation periods to discard (default 5).
- `λ_init`, `r_eval`, `abstol`, `reltol`: passed to `estimate_recovery_rate`
  when `method isa EnvelopeODE`.

# Returns
A `PolynomialSurface` with `flutter_frequency` set to `ω_c`.
"""
function forecast(
        dataset::TransientDataset,
        ω_c::Real;
        method::RecoveryRateMethod = EnvelopeODE(),
        r_grid::Union{Nothing, AbstractVector{<:Real}} = nothing,
        order::Int = 1,
        ridge::Real = 0.0,
        min_prominence::Real = 0.0,
        n_periods::Int = 5,
        λ_init::Union{Nothing, AbstractVector{<:Real}} = nothing,
        r_eval::Union{Nothing, AbstractVector{<:Real}} = nothing,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6
    )::PolynomialSurface

    isempty(dataset.samples) && error("forecast: dataset contains no samples")

    curves = RecoveryRateCurve[]

    for sample in dataset.samples
        # Phase fixing: extract local maxima, discard initial transient.
        t_pk, r_pk = find_peaks(sample.q, sample.t; min_prominence)
        t_pk, r_pk = discard_initial_transient(t_pk, r_pk, ω_c; n_periods)

        isempty(t_pk) && error(
            "forecast: no peaks remain after transient discard for μ=$(sample.μ). " *
                "Try reducing n_periods or increasing the signal length."
        )

        # Recovery-rate estimation — dispatch on method type.
        curve = if method isa EnvelopeODE
            kw = (abstol = abstol, reltol = reltol)
            if !isnothing(λ_init)
                kw = merge(kw, (λ_init = λ_init,))
            end
            if !isnothing(r_eval)
                kw = merge(kw, (r_eval = r_eval,))
            end
            estimate_recovery_rate(t_pk, r_pk, sample.μ, method; kw...)
        else
            estimate_recovery_rate(t_pk, r_pk, sample.μ, method)
        end

        push!(curves, curve)
    end

    # Build the common amplitude grid if not provided.
    r_out = if isnothing(r_grid)
        r_all = vcat([c.r for c in curves]...)
        range(minimum(r_all), maximum(r_all); length = 50) |> collect
    else
        collect(Float64, r_grid)
    end

    surface = fit_polynomial_surface(curves, r_out; order, ridge)

    # Attach flutter frequency.
    return PolynomialSurface(
        surface.r_grid, surface.coeffs, surface.param_mean,
        surface.param_scale, surface.order, surface.Np,
        surface.Nc, Float64(ω_c)
    )
end
