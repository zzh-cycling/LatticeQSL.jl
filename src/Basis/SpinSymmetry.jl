# ═══════════════════════════════════════════════════════════════
# SpinSymmetry: Sz conservation (U(1)) and spin inversion (Z₂)
# ═══════════════════════════════════════════════════════════════

"""
    SzConservation <: AbstractInternalSymmetry

Sz conservation symmetry: the total Sz = Σᵢ Sᵢᶻ is a good quantum number.
This is automatically handled by constructing `SpinBasis(nsites, nup)`.

# Fields
- `nup::Int` — number of up spins (total Sz = nup - N/2)
"""
struct SzConservation <: AbstractInternalSymmetry
    nup::Int
end

"""
    target_sz(sym::SzConservation, nsites::Int) -> Float64

Return the target Sz value.
"""
target_sz(sym::SzConservation, nsites::Int) = sym.nup - nsites / 2.0

"""
    SpinInversion <: AbstractInternalSymmetry

Spin inversion (Z₂) symmetry: P|s₁s₂⋯sN⟩ = |s̄₁s̄₂⋯s̄N⟩ where s̄ = 1-s.

This symmetry commutes with Sz=0 sector (it maps Sz=0 to itself).

# Fields
- `eigenvalue::Int` — +1 or -1, the eigenvalue of the spin inversion operator
"""
struct SpinInversion <: AbstractInternalSymmetry
    eigenvalue::Int
    function SpinInversion(eigenvalue::Int)
        @assert eigenvalue ∈ (-1, 1) "SpinInversion eigenvalue must be ±1"
        new(eigenvalue)
    end
end

"""
    invert_spins(state::BitStr{N, Int}) -> BitStr{N, Int}

Flip all spins: |↑↓↑⟩ → |↓↑↓⟩.
"""
function invert_spins(state::BitStr{N, Int}) where {N}
    T = BitStr{N, Int}
    mask = (1 << N) - 1
    return T(xor(state.buf, mask))
end
