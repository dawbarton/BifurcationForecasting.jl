# era.jl
#
# Eigensystem Realisation Algorithm (ERA) for modal basis extraction, plus
# biorthogonal modal projection and bandpass filtering.
#
# ERA constructs a minimal state-space model from multi-channel free-decay
# response data by SVD of the Hankel matrix.  For bifurcation forecasting the
# only output needed is the bifurcating mode's eigenvectors (the mode with the
# most positive real part), which are used to project subsequent transients onto
# a scalar modal coordinate.
#
# ERA is performed once, at the parameter sample closest to the flutter
# boundary.  The resulting ModalBasis is reused for all other samples; this is
# justified by the assumption that the centre-space eigenvectors vary slowly
# with the control parameter.
#
# References
# ──────────
#   Juang & Pappa (1985) J. Guid. Control Dyn. 8, 620–627  [original ERA]
#   Ghadami & Epureanu (2016) J. Comput. Nonlinear Dyn. 11, 061009, §3
#     https://doi.org/10.1115/1.4033920

# ── ERA ───────────────────────────────────────────────────────────────────────

"""
    era_modal_basis(t, Y; model_order, r_hankel, c_hankel) → ModalBasis

Extract the bifurcating-mode basis from multi-channel free-decay data using the
Eigensystem Realisation Algorithm.

# Arguments
- `t`: time vector (length N, uniform spacing assumed).
- `Y`: response matrix (n_outputs × N).

# Keyword arguments
- `model_order`: number of modes to retain in the truncated SVD (default 2;
  must be even as complex modes come in conjugate pairs).
- `r_hankel`: number of block rows in the Hankel matrix (default 100).
- `c_hankel`: number of block columns in the Hankel matrix (default 100).

# Returns
A `ModalBasis` whose eigenvalue has the most positive real part (least stable
mode, i.e. the bifurcating mode approaching flutter).

# Algorithm
1. Assemble Hankel matrices H(0) and H(1) from the output sequence.
2. SVD of H(0); truncate to `model_order` singular values.
3. Form discrete-time state-transition matrix S and output matrix C.
4. Diagonalise S; convert to continuous-time eigenvalues η = ln(s)/Δt.
5. Identify bifurcating mode (most positive Re(η)).
6. Form Φ = [Re(ϕ) Im(ϕ)] and Ψ = [Re(ψ) Im(ψ)]; pre-compute (ΨᵀΦ)⁻¹.
"""
function era_modal_basis(
        t::AbstractVector{<:Real},
        Y::AbstractMatrix{<:Real};
        model_order::Int = 2,
        r_hankel::Int = 100,
        c_hankel::Int = 100
    )::ModalBasis

    N = size(Y, 2)
    n_out = size(Y, 1)
    Δt = t[2] - t[1]

    @assert model_order % 2 == 0 "model_order must be even (conjugate pairs)"
    @assert r_hankel + c_hankel ≤ N "Hankel dimensions exceed data length"

    # Build block-Hankel matrices H(0) and H(1).
    # H(k) has block entry Y[:,j+k] at row j, column j.
    H0 = _hankel(Y, r_hankel, c_hankel, 0)
    H1 = _hankel(Y, r_hankel, c_hankel, 1)

    # Truncated SVD of H(0); retain first model_order singular triplets.
    F = svd(H0)
    P̄ = F.U[:, 1:model_order]        # (r_hankel*n_out × model_order)
    J̄ = F.V[:, 1:model_order]        # (c_hankel*n_out × model_order)
    Σ̄_sqrt = Diagonal(sqrt.(F.S[1:model_order]))
    Σ̄_isqrt = Diagonal(1.0 ./ sqrt.(F.S[1:model_order]))

    # Discrete-time state-transition and output matrices (Juang & Pappa 1985).
    S = Σ̄_isqrt * P̄' * H1 * J̄ * Σ̄_isqrt   # (model_order × model_order)
    C = (_block_first_row(n_out, r_hankel) * P̄ * Σ̄_sqrt)'  # (model_order × n_out) -> take transpose

    # Eigendecomposition of S.
    E = eigen(S)
    λ_disc = E.values    # discrete-time eigenvalues

    # Convert to continuous time: η = ln(s)/Δt.
    η = log.(complex.(λ_disc)) ./ Δt   # complex continuous-time eigenvalues

    # Identify bifurcating mode: most positive real part.
    # Eigenvalues come in conjugate pairs; pick the one with positive imaginary part.
    positive_idx = findall(imag.(η) .> 0)
    isempty(positive_idx) && error("ERA: no eigenvalue with positive imaginary part found")
    idx = positive_idx[argmax(real.(η[positive_idx]))]

    σ_c = η[idx]   # g_c + iω_c

    # Right eigenvector of S at index idx, then transform to original coordinates.
    ϕ_reduced = E.vectors[:, idx]   # (model_order,) complex
    # Ghadami 2016 eq. (6): physical output eigenvector = P̄ Σ̄^{1/2} ϕ_reduced.
    ϕ_phys = P̄ * Σ̄_sqrt * ϕ_reduced   # (r_hankel*n_out,) complex

    # Extract the first-block rows (length n_out) as the physical right mode shape.
    ϕ_c = ϕ_phys[1:n_out]   # (n_out,) complex

    # Left eigenvector: biorthogonal vector in the physical output space.
    # Assembling Φ_all (n_out × n_modes) with the first-block output eigenvector of
    # every positive-imaginary mode gives the output mode shape matrix.  The
    # biorthogonal left vectors satisfy ψ_k^H ϕ_j = δ_{k,j} and are the columns of
    # pinv(Φ_all^H).  For the square case (n_modes == n_out) this equals inv(Φ_all^H).
    # Biorthogonality in output space guarantees exact cancellation of each other mode
    # from the projected scalar coordinate q(t), whereas truncating the ERA state-space
    # left eigenvector to n_out rows does not preserve this property.
    n_modes = length(positive_idx)
    Φ_all = zeros(ComplexF64, n_out, n_modes)
    for (k, pidx) in enumerate(positive_idx)
        Φ_all[:, k] = (P̄ * Σ̄_sqrt * E.vectors[:, pidx])[1:n_out]
    end
    Ψ_all = pinv(Φ_all')          # (n_out × n_modes) biorthogonal left vectors
    bif_k = findfirst(==(idx), positive_idx)
    ψ_c = Ψ_all[:, bif_k]      # (n_out,) complex — biorthogonal to all other modes

    # Build real (N_phys × 2) matrices Φ and Ψ.
    Φ = [real(ϕ_c) imag(ϕ_c)]   # (n_out × 2)
    Ψ = [real(ψ_c) imag(ψ_c)]   # (n_out × 2)

    # Pre-compute (ΨᵀΦ)⁻¹ for use in modal projection.
    ΨᵀΦ_inv = inv(Ψ' * Φ)

    return ModalBasis(Φ, Ψ, ΨᵀΦ_inv, σ_c)
end

# ── Modal projection ──────────────────────────────────────────────────────────

"""
    modal_project(Y, basis::ModalBasis) → Vector{Float64}

Project a multi-channel response matrix Y (n_outputs × N) onto the scalar
modal coordinate of the bifurcating mode.

The biorthogonal projection is

    q(t) = [(ΨᵀΦ)⁻¹ Ψᵀ Δy(t)]₁

where the subscript 1 denotes the first component (real part of the modal
amplitude at phase θ = 0).

# Reference
Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, 103201, eq. (15).
"""
function modal_project(
        Y::AbstractMatrix{<:Real},
        basis::ModalBasis
    )::Vector{Float64}
    # Reconstruct complex vectors from the real column storage.
    ψ_c = basis.Ψ[:, 1] .+ im .* basis.Ψ[:, 2]   # (n_out,) complex
    ϕ_c = basis.Φ[:, 1] .+ im .* basis.Φ[:, 2]   # (n_out,) complex
    # Compute q(t) = Re(ψ_c^H Y(:,t) / ψ_c^H ϕ_c) directly in complex arithmetic.
    # ψ_c^H ϕ_c has magnitude ≈ 1 (well-conditioned), unlike the real 2×2 matrix
    # (Ψ^T Φ)^{-1} which becomes catastrophically ill-conditioned when ψ_c and ϕ_c
    # are nearly orthogonal in real space despite complex biorthogonality.
    scale = dot(conj.(ψ_c), ϕ_c)
    return [real(dot(conj.(ψ_c), Y[:, k]) / scale) for k in 1:size(Y, 2)]
end

# ── Bandpass filter ───────────────────────────────────────────────────────────

"""
    bandpass_filter(q, t, ω_c; bandwidth) → Vector{Float64}

Apply a zero-phase bandpass filter centred at ω_c (rad/time) to the scalar
signal q using FFT.

This is the simpler Option A from the algorithm documentation, appropriate
when other modes are well separated from ω_c.  For multi-mode systems with
closely-spaced frequencies, prefer ERA-based projection.

# Arguments
- `q`: scalar signal (length N).
- `t`: time vector (length N, uniform spacing assumed).
- `ω_c`: centre frequency in rad/time.
- `bandwidth`: half-bandwidth in rad/time (default 0.2 * ω_c, i.e. ±20%).
"""
function bandpass_filter(
        q::AbstractVector{<:Real},
        t::AbstractVector{<:Real},
        ω_c::Real;
        bandwidth::Real = 0.2 * ω_c
    )::Vector{Float64}

    N = length(q)
    Δt = t[2] - t[1]
    fs = 1.0 / Δt   # sample rate

    # FFT and corresponding frequencies (rad/time).
    Q = fft(q)
    freq = fftfreq(N, fs) .* 2π   # convert Hz → rad/time

    # Zero-phase: apply a rectangular window in the frequency domain.
    # Both positive and negative frequency bands are retained to preserve the
    # real-valued nature of the signal after IFFT.
    ω_lo = ω_c - bandwidth
    ω_hi = ω_c + bandwidth
    mask = (abs.(freq) .>= ω_lo) .& (abs.(freq) .<= ω_hi)

    Q_filtered = Q .* mask
    return real(ifft(Q_filtered))
end

# ── Private helpers ───────────────────────────────────────────────────────────

# Assemble a block-Hankel matrix with shift k.
# Each "block" is a single scalar column Y[:,j]; for n_out > 1 blocks are
# stacked vertically, giving a (r*n_out) × c matrix.
function _hankel(
        Y::AbstractMatrix{<:Real},
        r::Int, c::Int, k::Int
    )::Matrix{Float64}
    n_out = size(Y, 1)
    H = zeros(r * n_out, c)
    for col in 1:c, row in 1:r
        t_idx = row + col - 1 + k
        H[((row - 1) * n_out + 1):(row * n_out), col] = Y[:, t_idx]
    end
    return H
end

# Extract the identity matrix for the first n_out rows of a (r*n_out)-tall space.
function _block_first_row(n_out::Int, r::Int)::Matrix{Float64}
    E = zeros(n_out, r * n_out)
    E[:, 1:n_out] = I(n_out)
    return E
end
