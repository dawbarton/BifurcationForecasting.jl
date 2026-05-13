# signal.jl
#
# Signal processing utilities for the bifurcation forecasting pipeline.
#
# After modal projection (or bandpass filtering), the scalar signal q(t) is an
# oscillating, decaying waveform.  Phase fixing selects a consistent phase by
# extracting local maxima, giving a sparse (t_peak, r_peak) sequence that
# represents the amplitude envelope at fixed oscillation phase.
#
# An initial portion of the transient is discarded so that non-bifurcating
# modes have had time to decay and the data lies on the slow inertial manifold.
#
# Reference
# ─────────
#   Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009, §4.1–4.2
#     https://doi.org/10.1115/1.4033920

"""
    find_peaks(q, t; min_prominence) → (t_peaks, r_peaks)

Extract the times and amplitudes of local positive maxima of q(t).

A local maximum at index k satisfies q[k] > q[k-1] and q[k] > q[k+1].
Only positive maxima are retained (phase θ = 0⁺ in the oscillating signal).

# Arguments
- `q`: scalar signal (length N).
- `t`: corresponding time vector (length N).
- `min_prominence`: minimum peak height; peaks below this threshold are
  discarded (default 0.0, i.e. only positivity enforced).

# Returns
`(t_peaks, r_peaks)` — vectors of peak times and peak amplitudes.
"""
function find_peaks(
        q::AbstractVector{<:Real},
        t::AbstractVector{<:Real};
        min_prominence::Real = 0.0
    )

    t_peaks = Float64[]
    r_peaks = Float64[]

    for k in 2:(length(q) - 1)
        if q[k] > q[k - 1] && q[k] > q[k + 1] && q[k] > min_prominence
            push!(t_peaks, t[k])
            push!(r_peaks, q[k])
        end
    end

    return t_peaks, r_peaks
end

"""
    discard_initial_transient(t_peaks, r_peaks, ω_c; n_periods) → (t, r)

Remove peaks that fall within the first `n_periods` oscillation cycles.

This ensures that non-bifurcating modes (which decay faster than the
bifurcating mode) have become negligible before recovery-rate estimation
begins.  The oscillation period is T = 2π / ω_c.

# Arguments
- `t_peaks`, `r_peaks`: output of `find_peaks`.
- `ω_c`: flutter-mode angular frequency (rad/time).
- `n_periods`: number of initial oscillation periods to discard (default 5).
"""
function discard_initial_transient(
        t_peaks::AbstractVector{<:Real},
        r_peaks::AbstractVector{<:Real},
        ω_c::Real;
        n_periods::Int = 5
    )

    T_skip = n_periods * 2π / ω_c
    t0 = isempty(t_peaks) ? 0.0 : t_peaks[1]
    keep = t_peaks .>= t0 + T_skip

    return t_peaks[keep], r_peaks[keep]
end
