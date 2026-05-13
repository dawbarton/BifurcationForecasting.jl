# recovery_rate.jl
#
# Estimation of the recovery rate λ(r; μ) from a peak sequence.
#
# The recovery rate λ = ṙ/r measures how quickly the system returns to
# equilibrium from a given amplitude r.  It is negative in the pre-flutter
# regime and zero on the bifurcation diagram (flutter boundary or LCO branch).
#
# Two estimation strategies are provided:
#
#   FiniteDifference — direct log-ratio of consecutive peaks; requires many
#     peaks (> ~10) to give a reliable estimate across the amplitude range.
#
#   EnvelopeODE — fits a polynomial-in-r model for λ by integrating the
#     amplitude ODE ṙ = r (λ₀ + λ₁r + … + λₚrᵖ) and matching observed peaks
#     with Levenberg–Marquardt least squares.  Suitable when only a few peaks
#     (≥ 3) are available, as is common for low-frequency oscillatory systems.
#
# The EnvelopeODE approach uses NonlinearSolve.jl (LevenbergMarquardt solver)
# for the optimisation and OrdinaryDiffEq.jl for forward integration of the
# envelope ODE.
#
# References
# ──────────
#   Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009, §4.2–4.3
#     https://doi.org/10.1115/1.4033920
#   Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §2.2

# ── Finite-difference method ──────────────────────────────────────────────────

"""
    estimate_recovery_rate(t_peaks, r_peaks, μ, ::FiniteDifference) → RecoveryRateCurve

Estimate recovery rates from consecutive peak pairs via finite differences of
log(r):

    λ_k = (ln r_{k+1} − ln r_k) / (t_{k+1} − t_k)

Each consecutive pair yields one (r, λ) data point.  The finite-difference
ratio approximates λ at the midpoint time t_{k+1/2}, so the associated
amplitude is taken as the midpoint r_{k+1/2} = (r_k + r_{k+1})/2; this gives
O(Δt²) accuracy in the (r, λ) pairing rather than the O(Δt) that results from
associating with the left endpoint.
"""
function estimate_recovery_rate(
        t_peaks::AbstractVector{<:Real},
        r_peaks::AbstractVector{<:Real},
        μ::AbstractVector{<:Real},
        ::FiniteDifference
    )::RecoveryRateCurve

    n = length(t_peaks)
    n < 2 && error("FiniteDifference requires at least 2 peaks; got $n")

    r_out = zeros(n - 1)
    λ_out = zeros(n - 1)

    for k in 1:(n - 1)
        Δt = t_peaks[k + 1] - t_peaks[k]
        λ_out[k] = (log(r_peaks[k + 1]) - log(r_peaks[k])) / Δt
        # Midpoint amplitude: the FD ratio approximates λ at the midpoint time
        # t_{k+1/2}, so associating it with (r_k + r_{k+1})/2 gives O(Δt²)
        # rather than O(Δt) accuracy in the (r, λ) pairing.
        r_out[k] = (r_peaks[k] + r_peaks[k + 1]) / 2
    end

    return RecoveryRateCurve(collect(μ), r_out, λ_out, FiniteDifference())
end

# ── Envelope-ODE method ───────────────────────────────────────────────────────

"""
    estimate_recovery_rate(t_peaks, r_peaks, μ, method::EnvelopeODE;
                           λ_init, r_eval) → RecoveryRateCurve

Fit the amplitude-envelope ODE

    ṙ = r (λ₀ + λ₁r + λ₂r² + … + λₚrᵖ)                     (*)

to observed peak amplitudes (t_peaks, r_peaks) using Levenberg–Marquardt
nonlinear least squares (NonlinearSolve.jl).

Once the coefficients [λ₀, …, λₚ] are found, the recovery rate is evaluated
at a dense set of amplitude values spanning [min(r_peaks), max(r_peaks)].

# Keyword arguments
- `λ_init`: initial guess for coefficients (length p+1); defaults to
  `[g_c, 0, …, 0]` if the linear damping g_c = real(σ_c) is available, or
  `zeros(p+1)` otherwise.
- `r_eval`: amplitude grid for output evaluation (default: 50 points, linear).
- `ode_dt`: fixed step for envelope ODE integration (default 0.1 × mean peak
  spacing; use a small value relative to the peak separation).
- `abstol`, `reltol`: ODE solver tolerances (defaults 1e-8, 1e-6).
"""
function estimate_recovery_rate(
        t_peaks::AbstractVector{<:Real},
        r_peaks::AbstractVector{<:Real},
        μ::AbstractVector{<:Real},
        method::EnvelopeODE;
        λ_init::Union{Nothing, AbstractVector{<:Real}} = nothing,
        r_eval::Union{Nothing, AbstractVector{<:Real}} = nothing,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6
    )::RecoveryRateCurve

    p = method.poly_order
    n_peaks = length(t_peaks)
    n_peaks < 3 && error("EnvelopeODE requires at least 3 peaks; got $n_peaks")

    # Default initial guess: log-slope of the peak sequence for λ₀, zeros for
    # higher-order terms.  For a decaying signal this gives λ₀ ≈ g_c, which is
    # the correct sign and magnitude and avoids the spurious local minimum at
    # λ₀ > 0 that the zero initial guess reaches for finite-difference Jacobians.
    λ0 = if isnothing(λ_init)
        λ0_init = (log(r_peaks[end]) - log(r_peaks[1])) / (t_peaks[end] - t_peaks[1])
        [λ0_init; zeros(p)]
    else
        collect(Float64, λ_init)
    end
    length(λ0) == p + 1 || error("λ_init must have length poly_order+1 = $(p + 1)")

    # Residual function for NonlinearSolve: returns vector of r_model(t_k) - r_obs_k.
    # Closure captures ODE solver settings and observed data.
    function residual!(res, λ_vec, _)
        r_model = _integrate_envelope(λ_vec, t_peaks, r_peaks[1]; abstol, reltol)
        return @. res = r_model - r_peaks
    end

    prob = NonlinearLeastSquaresProblem(
        NonlinearFunction(residual!; resid_prototype = zeros(n_peaks)),
        λ0
    )
    sol = solve(prob, LevenbergMarquardt(; autodiff = AutoFiniteDiff()); maxiters = 500)

    λ_fit = sol.u

    # Evaluate the fitted λ(r) = λ₀ + λ₁r + … + λₚrᵖ on a dense amplitude grid.
    r_min = minimum(r_peaks)
    r_max = maximum(r_peaks)
    r_out = isnothing(r_eval) ? range(r_min, r_max; length = 50) |> collect :
        collect(Float64, r_eval)

    λ_out = [_eval_poly(λ_fit, r) for r in r_out]

    return RecoveryRateCurve(collect(μ), r_out, λ_out, method)
end

# ── Private helpers ───────────────────────────────────────────────────────────

# Evaluate the polynomial λ₀ + λ₁r + … + λₚrᵖ at scalar r (Horner's method).
function _eval_poly(coeffs::AbstractVector{<:Real}, r::Real)::Float64
    val = coeffs[end]
    for k in (length(coeffs) - 1):-1:1
        val = val * r + coeffs[k]
    end
    return val
end

# Integrate the envelope ODE ṙ = r * poly(r) from r0, evaluating at t_eval.
# Uses OrdinaryDiffEq with Tsit5 (explicit, RK-based).
function _integrate_envelope(
        λ_vec::AbstractVector{<:Real},
        t_eval::AbstractVector{<:Real},
        r0::Real;
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6
    )::Vector{Float64}

    # Envelope ODE: ṙ = r (λ₀ + λ₁r + … + λₚrᵖ)
    function envelope_ode!(du, u, p, t)
        r = u[1]
        return du[1] = r * _eval_poly(p, r)
    end

    tspan = (Float64(t_eval[1]), Float64(t_eval[end]))
    prob = ODEProblem(envelope_ode!, [Float64(r0)], tspan, λ_vec)
    sol = solve(
        prob, Tsit5();
        saveat = collect(Float64, t_eval),
        abstol = abstol,
        reltol = reltol
    )

    # Extract first state component at each requested time; return NaN if solve failed.
    if sol.retcode == ReturnCode.Success
        return [sol(t)[1] for t in t_eval]
    else
        return fill(NaN, length(t_eval))
    end
end
