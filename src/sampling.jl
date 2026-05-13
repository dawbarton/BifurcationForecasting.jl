# sampling.jl
#
# Utility functions for generating parameter grids and amplitude grids.
#
# Latin Hypercube Sampling (LHS) distributes Ns samples across an Np-dimensional
# parameter box such that each marginal dimension is covered uniformly, while
# random permutations break the aligned-grid regularity.  This gives better
# space-filling properties than a uniform grid with the same number of points.
#
# References
# ──────────
#   McKay, Beckman & Conover (1979) Technometrics 21, 239–245  [LHS original]
#   Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, §3.1
#     https://doi.org/10.1016/j.jfluidstructs.2020.103201

"""
    lhs_samples(Ns, lb, ub; rng) → Vector{Vector{Float64}}

Generate `Ns` parameter samples by Latin Hypercube Sampling in the box
defined by lower bounds `lb` and upper bounds `ub` (both length-Np vectors).

Each dimension is divided into `Ns` equal strata.  One sample is drawn
uniformly from each stratum, and the samples are randomly permuted
independently per dimension.  This ensures that the marginal distribution
in every parameter is uniform while avoiding the lattice structure of a
regular grid.

# Arguments
- `Ns`: number of samples.
- `lb`: lower bounds (length Np).
- `ub`: upper bounds (length Np).

# Keyword arguments
- `rng`: random number generator (default `Random.default_rng()`).

# Returns
A vector of `Ns` parameter vectors, each of length `Np`.
"""
function lhs_samples(
        Ns::Int,
        lb::AbstractVector{<:Real},
        ub::AbstractVector{<:Real};
        rng = Random.default_rng()
    )::Vector{Vector{Float64}}

    Np = length(lb)
    length(ub) == Np || error("lb and ub must have the same length")
    all(ub .> lb) || error("ub must be strictly greater than lb in every dimension")

    # Stratum width per dimension.
    h = (collect(Float64, ub) .- collect(Float64, lb)) ./ Ns

    # (Ns × Np) matrix: row l is the l-th sample, column i is parameter i.
    X = zeros(Ns, Np)
    for i in 1:Np
        perm = randperm(rng, Ns)
        for l in 1:Ns
            # Draw uniformly from stratum perm[l] of dimension i.
            X[l, i] = lb[i] + h[i] * (perm[l] - 1 + rand(rng))
        end
    end

    return [X[l, :] for l in 1:Ns]
end

"""
    amplitude_grid(Nr, r_max; r_min, spacing) → Vector{Float64}

Generate an amplitude grid of `Nr` values in the range `[r_min, r_max]`.

# Arguments
- `Nr`: number of grid points.
- `r_max`: maximum amplitude (sets the upper bound; must exceed `r_min`).

# Keyword arguments
- `r_min`: minimum amplitude (default `r_max / Nr`; a small positive value
  to avoid r = 0 singularities in recovery-rate evaluation).
- `spacing`: `:linear` (default) or `:sqrt` for denser coverage near r = 0.
  Square-root spacing improves flutter-boundary resolution.

# Returns
Sorted vector of amplitude values in ascending order.

# Example
```julia
r_grid = amplitude_grid(40, 0.20; r_min = 0.001, spacing = :sqrt)
```
"""
function amplitude_grid(
        Nr::Int,
        r_max::Real;
        r_min::Real = r_max / Nr,
        spacing::Symbol = :linear
    )::Vector{Float64}

    Nr > 0   || error("Nr must be positive")
    r_max > 0 || error("r_max must be positive")
    r_min > 0 || error("r_min must be positive; use the default or set a small value > 0")
    r_max > r_min || error("r_max must exceed r_min")

    if spacing == :linear
        return collect(range(Float64(r_min), Float64(r_max); length = Nr))
    elseif spacing == :sqrt
        # Map linearly in √r so that small amplitudes are more densely sampled.
        grid_sqrt = range(sqrt(Float64(r_min)), sqrt(Float64(r_max)); length = Nr)
        return collect(x^2 for x in grid_sqrt)
    else
        error("spacing must be :linear or :sqrt; got :$spacing")
    end
end
