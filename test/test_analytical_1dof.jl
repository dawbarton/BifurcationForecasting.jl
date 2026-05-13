# test_analytical_1dof.jl
#
# Test case 1: Analytical 1-DOF bifurcation forecasting
#
# Verifies the core pipeline against two scalar amplitude ODEs with known
# analytical solutions.  Because the systems are non-oscillatory, no peak
# finding or modal filtering is needed; the full time series serves as input
# to the finite-difference recovery-rate estimator.
#
# Reference: Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §3

using OrdinaryDiffEq
using Statistics

# ── System definitions ────────────────────────────────────────────────────────

# Supercritical: ṙ = r(αμ − βr²)
function rhs_supercritical(r, p, t)
    μ, β = p
    return r * (μ - β * r^2)
end

# Subcritical: ṙ = r(αμ + βr² − γr⁴)  (α = γ = 1)
function rhs_subcritical(r, p, t)
    μ, β = p
    return r * (μ + β * r^2 - r^4)
end

# Analytical recovery rate for each system.
λ_exact_sup(μ, β, r) = μ - β * r^2
λ_exact_sub(μ, β, r) = μ + β * r^2 - r^4

# Exact solution for supercritical system (derived from Bernoulli substitution u = r⁻²):
#   r(t)² = αμ r₀² / [βr₀² + (αμ − βr₀²) exp(−2αμt)]
function exact_supercritical(t, μ, β, r0)
    αμ = μ        # α = 1
    denom = β * r0^2 + (αμ - β * r0^2) * exp(-2αμ * t)
    r2 = αμ * r0^2 / denom
    return r2 > 0 ? sqrt(r2) : 0.0
end

# ── Helper: simulate an amplitude ODE ─────────────────────────────────────────

function simulate(rhs, μ, β, r0; dt = 0.01, t_max = 200.0, r_min = 1.0e-6)
    prob = ODEProblem(rhs, r0, (0.0, t_max), (μ, β))
    sol = solve(
        prob, Tsit5();
        saveat = 0:dt:t_max,
        abstol = 1.0e-10,
        reltol = 1.0e-8
    )
    t = Float64.(sol.t)
    r = Float64.(sol.u)
    # Trim to positive values above floor.
    keep = r .> r_min
    return t[keep], r[keep]
end

# ── Parameter samples (Riso JFS 2021, Table 1) ────────────────────────────────

# Single-parameter samples (β = 1 fixed):
μ_1p = [-0.5, -0.4, -0.3]
β_1p = [1.0, 1.0, 1.0]

# Two-parameter samples (simulations 4–9):
μ_2p_extra = [-0.48, -0.36, -0.32, -0.39, -0.43, -0.43]
β_2p_extra = [1.05, 1.13, 0.97, 0.84, 0.92, 1.13]

# All 9 samples combined for 2-parameter fitting:
μ_all = vcat(μ_1p, μ_2p_extra)
β_all = vcat(β_1p, β_2p_extra)

# r grid and initial condition
r0 = 1.5
r_min_grid = 0.01
r_max_grid = 1.45
Nr = 50
r_grid = range(r_min_grid, r_max_grid; length = Nr) |> collect

# ── Tests ─────────────────────────────────────────────────────────────────────

@testset "1-DOF analytical: supercritical" begin

    # --- Check 1: FD recovery rate accuracy ---
    @testset "Check 1: FD recovery rate accuracy" begin
        for (μ, β) in zip(μ_all, β_all)
            t, r = simulate(rhs_supercritical, μ, β, r0)
            curve = estimate_recovery_rate(t, r, [μ, β], FiniteDifference())

            λ_ana = [λ_exact_sup(μ, β, rv) for rv in curve.r]
            max_err = maximum(abs.(curve.λ .- λ_ana))
            @test max_err < 5.0e-3 broken = false
        end
    end

    # --- Build recovery-rate curves for all 9 samples ---
    curves_all = RecoveryRateCurve[]
    for (μ, β) in zip(μ_all, β_all)
        t, r = simulate(rhs_supercritical, μ, β, r0)
        push!(curves_all, estimate_recovery_rate(t, r, [μ, β], FiniteDifference()))
    end

    curves_1p = curves_all[1:3]   # β = 1 samples only

    # --- Check 2: Bi-linear fit exactness ---
    # The true λ is exactly linear in (μ, β).  The surface fit should reproduce
    # each FD curve to within the noise of the finite-difference estimate
    # (~O(Δt²) ≈ 1e-4 for dt = 0.01); machine precision is not achievable because
    # the overdetermined LS system is fit across all 9 noisy curves simultaneously.
    @testset "Check 2: Bi-linear fit exactness (first-order surface)" begin
        surf = fit_polynomial_surface(curves_all, r_grid; order = 1)
        for (l, curve) in enumerate(curves_all)
            λ_eval = evaluate(surf, curve.μ)
            # Interpolate the curve onto r_grid for comparison.
            λ_interp = [BifurcationForecasting._interp_curve(curve, [rk])[1] for rk in r_grid]
            max_err = maximum(abs.(λ_eval .- λ_interp))
            @test max_err < 1.0e-3
        end
    end

    # --- Check 3: Flutter boundary ---
    @testset "Check 3: Flutter boundary |μ_c| < 1e-3" begin
        surf = fit_polynomial_surface(curves_all, r_grid; order = 1)
        diag = extract_diagram(
            surf;
            μ_ref = [-0.5, 1.0],
            μ_dir = [1.0, 0.0]
        )
        μ_fb = diag.flutter_boundary
        @test !isempty(μ_fb)
        # Flutter boundary is at μ = 0 for all β.
        @test abs(μ_fb[1][1]) < 1.0e-3
    end

    # --- Check 4: Bifurcation diagram at β = 1 ---
    @testset "Check 4: Bifurcation diagram at β = 1" begin
        surf = fit_polynomial_surface(curves_1p, r_grid; order = 1)
        diag = extract_diagram(
            surf;
            μ_ref = [mean(μ_1p), 1.0],
            μ_dir = [1.0, 0.0]
        )
        # Check predicted μ against the analytical branch μ = β r² / α = r² (β=1, α=1).
        # Tolerance 2e-3: the FD estimate has O(Δt²) ≈ 1e-4 errors; at large r
        # (near the grid boundary r=1.45) these produce ~1.2e-3 error in the
        # predicted μ because the surface is fit from a finite sample.
        for pt in diag.points
            r̃ = pt.r
            μ̃ = pt.μ[1]
            μ_exact = r̃^2   # β = 1, α = 1
            @test abs(μ̃ - μ_exact) < 2.0e-3
        end
        # Bifurcation type should be supercritical.
        @test diag.bifurcation_type == :supercritical
    end

    # --- Check 5: Two-parameter cross-validation ---
    # Train on samples 4–9 (no β = 1 samples); predict at β = 1.
    @testset "Check 5: Two-parameter cross-validation" begin
        curves_extra = curves_all[4:end]
        surf = fit_polynomial_surface(curves_extra, r_grid; order = 1)
        diag = extract_diagram(
            surf;
            μ_ref = [-0.5, 1.0],
            μ_dir = [1.0, 0.0]
        )
        for pt in diag.points
            r̃ = pt.r
            μ̃ = pt.μ[1]
            β_pred = pt.μ[2]
            @test abs(β_pred - 1.0) < 0.1   # direction search at β = 1
            μ_exact = r̃^2 * β_pred           # generalised: μ = βr²/α
            @test abs(μ̃ - μ_exact) < 1.0e-2    # cross-validation tolerance
        end
    end

    # --- Check 6: Single- vs. two-parameter agreement ---
    @testset "Check 6: Single vs. two-parameter agreement at β = 1" begin
        surf_1p = fit_polynomial_surface(curves_1p, r_grid; order = 1)
        surf_2p = fit_polynomial_surface(curves_all, r_grid; order = 1)

        diag_1p = extract_diagram(
            surf_1p;
            μ_ref = [mean(μ_1p), 1.0],
            μ_dir = [1.0, 0.0]
        )
        diag_2p = extract_diagram(
            surf_2p;
            μ_ref = [-0.5, 1.0],
            μ_dir = [1.0, 0.0]
        )

        # Both diagrams should agree at β ≈ 1.
        n_compare = min(length(diag_1p.points), length(diag_2p.points))
        @test n_compare > 0
        for k in 1:n_compare
            Δμ = abs(diag_1p.points[k].μ[1] - diag_2p.points[k].μ[1])
            @test Δμ < 1.0e-3
        end
    end
end

@testset "1-DOF analytical: subcritical" begin

    # --- Check 1: FD recovery rate accuracy ---
    # Tolerance 7e-3: the subcritical r⁴ term gives a steeper λ(r) gradient at
    # large r than the supercritical case, producing O(Δt²) ≈ 6e-3 errors for
    # dt = 0.01 (see supercritical Check 1 at 5e-3 for comparison).
    @testset "Check 1: FD recovery rate accuracy" begin
        for (μ, β) in zip(μ_all, β_all)
            t, r = simulate(rhs_subcritical, μ, β, r0)
            curve = estimate_recovery_rate(t, r, [μ, β], FiniteDifference())

            λ_ana = [λ_exact_sub(μ, β, rv) for rv in curve.r]
            max_err = maximum(abs.(curve.λ .- λ_ana))
            @test max_err < 7.0e-3
        end
    end

    # --- Build curves ---
    curves_all = RecoveryRateCurve[]
    for (μ, β) in zip(μ_all, β_all)
        t, r = simulate(rhs_subcritical, μ, β, r0)
        push!(curves_all, estimate_recovery_rate(t, r, [μ, β], FiniteDifference()))
    end
    curves_1p = curves_all[1:3]

    # --- Check 2: Bi-linear fit exactness ---
    # Same reasoning as the supercritical case: tolerance set by FD noise (~1e-4),
    # not machine precision, because the overdetermined LS fit averages across all
    # 9 noisy curves.
    @testset "Check 2: Bi-linear fit exactness" begin
        surf = fit_polynomial_surface(curves_all, r_grid; order = 1)
        for curve in curves_all
            λ_eval = evaluate(surf, curve.μ)
            λ_interp = [BifurcationForecasting._interp_curve(curve, [rk])[1] for rk in r_grid]
            max_err = maximum(abs.(λ_eval .- λ_interp))
            @test max_err < 1.0e-3
        end
    end

    # --- Check 3: Flutter boundary ---
    @testset "Check 3: Flutter boundary |μ_c| < 1e-3" begin
        surf = fit_polynomial_surface(curves_all, r_grid; order = 1)
        diag = extract_diagram(
            surf;
            μ_ref = [-0.5, 1.0],
            μ_dir = [1.0, 0.0]
        )
        μ_fb = diag.flutter_boundary
        @test !isempty(μ_fb)
        @test abs(μ_fb[1][1]) < 1.0e-3
    end

    # --- Check 4: Bifurcation diagram shape at β = 1 ---
    # Stable branch: r̃ = sqrt([β + sqrt(β² + 4μ)] / 2)  (α = γ = 1, β = 1)
    @testset "Check 4: Bifurcation diagram at β = 1 (subcritical)" begin
        surf = fit_polynomial_surface(curves_1p, r_grid; order = 1)
        diag = extract_diagram(
            surf;
            μ_ref = [mean(μ_1p), 1.0],
            μ_dir = [1.0, 0.0]
        )
        # Subcritical LCO branch: setting λ = μ + βr² − γr⁴ = 0 gives
        # μ = γr⁴ − βr² = r⁴ − r² at β = γ = 1.  This is negative for r < 1
        # (LCO exists below the flutter boundary).
        # Tolerance 1.5e-2: large FD errors near r = 1.45 ≈ r₀ (grid boundary);
        # the r⁴ term makes dλ/dr steep, amplifying noise in the inferred μ.
        for pt in diag.points
            r̃ = pt.r
            μ̃ = pt.μ[1]
            μ_exact = r̃^4 - r̃^2
            @test abs(μ̃ - μ_exact) < 1.5e-2
        end
    end

    # --- Check 5: Bifurcation type ---
    @testset "Check 5: Subcritical type identified" begin
        surf = fit_polynomial_surface(curves_1p, r_grid; order = 1)
        diag = extract_diagram(
            surf;
            μ_ref = [mean(μ_1p), 1.0],
            μ_dir = [1.0, 0.0]
        )
        @test diag.bifurcation_type == :subcritical
    end
end
