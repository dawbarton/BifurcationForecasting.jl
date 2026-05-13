# state_velocity.jl
#
# Model-based recovery-rate estimation via the state velocity method.
#
# Instead of integrating free-decay transients, recovery rates are computed
# directly by evaluating the full nonlinear vector field f(μ, y) at states
# perturbed along the bifurcating mode eigenvector.  This requires a model
# (state velocity function and Jacobian), but avoids costly time integration.
#
# For each parameter value μᵢ and amplitude rₖ:
#
#   1. Find equilibrium yₑ satisfying f(μᵢ, yₑ) = 0.
#   2. Compute Jacobian Aᵢ = ∂f/∂y|_{μᵢ, yₑ}; identify bifurcating mode.
#   3. Form Φᵢ = [Re(ϕ_c) Im(ϕ_c)], Ψᵢ = [Re(ψ_c) Im(ψ_c)].
#   4. Perturbed state: y₀(rₖ) = yₑ + rₖ Re(ϕ_c).
#   5. Nonlinear part: f_nl = f(μᵢ, y₀) − Aᵢ (y₀ − yₑ).
#   6. Recovery rate:
#        λᵢₖ = gᵢ + [(Ψᵢᵀ Φᵢ)⁻¹ Ψᵢᵀ f_nl(y₀(rₖ))]₁ / rₖ
#
# Public API
# ──────────
#   forecast_state_velocity(f, J, μ_samples, r_grid; equilibria, order, ridge)
#       → PolynomialSurface
#
# Reference
# ─────────
#   Riso, Cesnik & Epureanu (2022) AIAA J. 60, 5401–5413, §2, eq. (20)
#     https://doi.org/10.2514/1.J061860

"""
    forecast_state_velocity(f, J, μ_samples, r_grid;
                            equilibria, flutter_frequency, order, ridge)
        → PolynomialSurface

Estimate the recovery-rate surface using the state velocity method (Riso 2022).

No transient time integration is required; recovery rates are obtained by
evaluating the nonlinear vector field at states perturbed along the bifurcating
mode eigenvector and projecting back onto the reduced modal space.

# Arguments
- `f`: state velocity function with signature `f(μ, y) → ẏ` (a vector of the
  same length as `y`).
- `J`: Jacobian function `J(μ, y_e) → A` returning the state matrix at
  equilibrium (`n × n` matrix).
- `μ_samples`: vector of control-parameter vectors (length Ns), each of length Np.
- `r_grid`: amplitude values at which to evaluate λ (length Nr).

# Keyword arguments
- `equilibria`: pre-computed equilibrium states (vector of `n`-vectors, one per
  parameter sample).  Default: zero vector (appropriate for systems whose
  equilibrium is at the origin for all μ).
- `flutter_frequency`: ω_c to store in the returned surface (default: taken from
  the bifurcating mode at the first sample).
- `order`: Taylor expansion order for the polynomial surface fit (1 or 2; default 1).
- `ridge`: Tikhonov regularisation for the surface fit (default 0.0).

# Returns
A `PolynomialSurface` identical in structure to the output of `forecast`.

# Reference
Riso, Cesnik & Epureanu (2022) AIAA J. 60, 5401–5413, §2 (eq. 20/25).
"""
function forecast_state_velocity(
        f,
        J,
        μ_samples::AbstractVector,
        r_grid::AbstractVector{<:Real};
        equilibria::Union{Nothing, AbstractVector} = nothing,
        n_states::Union{Nothing, Int} = nothing,
        flutter_frequency::Union{Nothing, Real} = nothing,
        order::Int = 1,
        ridge::Real = 0.0
    )::PolynomialSurface

    Ns = length(μ_samples)
    Ns < 1 && error("forecast_state_velocity: need at least one parameter sample")

    # Determine state dimension from equilibria or n_states.
    n = if !isnothing(equilibria)
        length(equilibria[1])
    elseif !isnothing(n_states)
        n_states
    else
        error("Provide either `equilibria` or `n_states` (for origin-equilibrium systems)")
    end

    r_grid_f = collect(Float64, r_grid)
    Nr = length(r_grid_f)

    curves = Vector{RecoveryRateCurve}(undef, Ns)

    ω_c_est = isnothing(flutter_frequency) ? NaN : Float64(flutter_frequency)

    for l in 1:Ns
        μ_l = collect(Float64, μ_samples[l])
        y_e = isnothing(equilibria) ? zeros(n) : collect(Float64, equilibria[l])

        # State matrix and bifurcating mode eigenvectors.
        A_l = J(μ_l, y_e)
        Φ_l, Ψ_l, ΨᵀΦ_inv_l, σ_c = _sv_bifurcating_mode(A_l)
        g_c = real(σ_c)

        # Store flutter frequency from first sample if not provided.
        if isnan(ω_c_est) && l == 1
            ω_c_est = abs(imag(σ_c))
        end

        # ϕ_re = Re(ϕ_c) = first column of Φ_l (physical space, n-vector).
        ϕ_re = Φ_l[:, 1]

        # Evaluate recovery rate at each amplitude on the grid.
        λ_l = zeros(Nr)
        for k in 1:Nr
            rₖ = r_grid_f[k]
            y0 = y_e .+ rₖ .* ϕ_re

            # Full state velocity and its nonlinear part (subtract linear contribution).
            f_full = f(μ_l, y0)
            f_nl = f_full .- A_l * (y0 .- y_e)

            # Biorthogonal projection: (Ψᵀ Φ)⁻¹ Ψᵀ f_nl  → (2,) vector.
            # The first component gives the nonlinear contribution to ṙ.
            proj = ΨᵀΦ_inv_l * (Ψ_l' * f_nl)

            # Recovery rate = linear part (g_c) + nonlinear part / r.
            λ_l[k] = g_c + proj[1] / rₖ
        end

        curves[l] = RecoveryRateCurve(μ_l, r_grid_f, λ_l, FiniteDifference())
    end

    surface = fit_polynomial_surface(curves, r_grid_f; order = order, ridge = ridge)

    # Attach flutter frequency (PolynomialSurface is immutable; rebuild).
    return PolynomialSurface(
        surface.r_grid, surface.coeffs, surface.param_mean,
        surface.param_scale, surface.order, surface.Np,
        surface.Nc, ω_c_est
    )
end

# ── Private helpers ───────────────────────────────────────────────────────────

# Extract bifurcating mode from an n × n state matrix A.
# Returns (Φ, Ψ, (ΨᵀΦ)⁻¹, σ_c) where σ_c is the eigenvalue with the most
# positive real part and positive imaginary part (least stable mode).
function _sv_bifurcating_mode(A::AbstractMatrix{<:Real})
    E = eigen(A)
    η = E.values   # complex eigenvalues

    # Identify the bifurcating mode: most positive Re(η), with Im(η) > 0.
    pos_idx = findall(imag.(η) .> 0)
    isempty(pos_idx) && error("_sv_bifurcating_mode: no eigenvalue with positive imaginary part")
    idx = pos_idx[argmax(real.(η[pos_idx]))]

    σ_c = η[idx]
    ϕ_c = E.vectors[:, idx]            # right eigenvector (complex, n-vector)
    ψ_c_raw = conj(inv(E.vectors)'[:, idx]) # left eigenvector

    # Real (n × 2) matrices of real and imaginary parts.
    Φ = [real(ϕ_c) imag(ϕ_c)]
    Ψ = [real(ψ_c_raw) imag(ψ_c_raw)]

    ΨᵀΦ_inv = inv(Ψ' * Φ)

    return Φ, Ψ, ΨᵀΦ_inv, σ_c
end
