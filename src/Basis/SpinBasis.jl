# ═══════════════════════════════════════════════════════════════
# SpinBasis: spin-1/2 basis using BitBasis.DitStr
#
# Each computational basis state |s₁ s₂ ⋯ sN⟩ is represented as a
# DitStr{2, N, Int}, where sᵢ ∈ {0, 1} (spin down/up).
#
# The full Hilbert space has dimension 2^N.
# With Sz conservation, we keep only states with a fixed number of up-spins.
# ═══════════════════════════════════════════════════════════════

using BitBasis: DitStr, BitStr

"""
    SpinBasis{N}

A basis for `N` spin-1/2 sites. States are stored as a **sorted** `Vector{BitStr{N, Int}}`
for efficient lookup via `searchsorted`.

# Fields
- `nsites::Int` — number of sites N
- `states::Vector{BitStr{N, Int}}` — sorted list of basis states
- `dim::Int` — dimension of the basis (= length(states))

# Type alias
`BitStr{N, Int} = DitStr{2, N, Int}` — a `DitStr` with dit = 2.
"""
struct SpinBasis{N}
    nsites::Int
    states::Vector{BitStr{N, Int}}
    dim::Int
end

"""
    SpinBasis(nsites::Int)

Construct the full (no symmetry) spin-1/2 basis for `nsites` sites.
Dimension = 2^nsites.
"""
function SpinBasis(nsites::Int)
    N = nsites
    T = BitStr{N, Int}
    states = [T(i) for i in 0:(2^N - 1)]
    return SpinBasis{N}(N, states, length(states))
end

"""
    SpinBasis(nsites::Int, nup::Int)

Construct the Sz-conserving spin-1/2 basis for `nsites` sites
with exactly `nup` up-spins (Sz = nup - nsites/2).
Dimension = C(nsites, nup).
"""
function SpinBasis(nsites::Int, nup::Int)
    N = nsites
    @assert 0 <= nup <= N "nup must be in [0, $N], got $nup"
    T = BitStr{N, Int}
    states = T[]
    for i in 0:(2^N - 1)
        count_ones(i) == nup && push!(states, T(i))
    end
    sort!(states)
    return SpinBasis{N}(N, states, length(states))
end

"""
    state_index(basis::SpinBasis{N}, state::BitStr{N, Int}) -> Int

Find the index of `state` in the basis. Returns 0 if not found.
Uses binary search (O(log dim)).
"""
function state_index(basis::SpinBasis{N}, state::BitStr{N, Int}) where {N}
    idx = searchsortedfirst(basis.states, state)
    if idx <= basis.dim && basis.states[idx] == state
        return idx
    end
    return 0
end

"""
    Base.getindex(basis::SpinBasis, i::Int)

Get the i-th basis state.
"""
Base.getindex(basis::SpinBasis, i::Int) = basis.states[i]
Base.length(basis::SpinBasis) = basis.dim
Base.iterate(basis::SpinBasis, args...) = iterate(basis.states, args...)

# ─────────── Bit manipulation helpers ───────────

"""
    spin_up(state::BitStr, site::Int) -> Bool

Check if spin at `site` (1-based) is up (= 1).
"""
spin_up(state::BitStr{N, Int}, site::Int) where {N} = ((state.buf >> (site - 1)) & 0x1) == 0x1

"""
    spin_value(state::BitStr, site::Int) -> Int

Return the spin value at `site` (1-based): +1 for up, -1 for down.
(For Sz eigenvalue: multiply by 1/2)
"""
spin_value(state::BitStr{N, Int}, site::Int) where {N} = 2 * Int((state.buf >> (site - 1)) & 0x1) - 1

"""
    flip_spin(state::BitStr{N, Int}, site::Int) -> BitStr{N, Int}

Flip the spin at `site` (1-based).
"""
function flip_spin(state::BitStr{N, Int}, site::Int) where {N}
    T = BitStr{N, Int}
    return T(xor(state.buf, 1 << (site - 1)))
end

"""
    total_sz(state::BitStr{N, Int}) -> Float64

Compute the total Sz = (number of up spins) - N/2.
"""
total_sz(state::BitStr{N, Int}) where {N} = count_ones(state.buf) - N / 2.0

"""
    num_up(state::BitStr{N, Int}) -> Int

Count the number of up spins.
"""
num_up(state::BitStr{N, Int}) where {N} = count_ones(state.buf)

"""
    permute_bits(state::BitStr{N, Int}, perm::Vector{Int}) -> BitStr{N, Int}

Apply a site permutation to a basis state.
`perm[i] = j` means: the spin at site `i` in the new state comes from site `j` in the old state.
"""
function permute_bits(state::BitStr{N, Int}, perm::Vector{Int}) where {N}
    T = BitStr{N, Int}
    result = 0
    for i in 1:N
        if ((state.buf >> (perm[i] - 1)) & 0x1) == 0x1
            result |= (1 << (i - 1))
        end
    end
    return T(result)
end
