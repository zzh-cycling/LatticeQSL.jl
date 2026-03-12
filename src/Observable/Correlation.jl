# ═══════════════════════════════════════════════════════════════
# Spin-spin correlations
# ═══════════════════════════════════════════════════════════════

using LinearAlgebra

"""
    szsz_correlation(basis::SpinBasis{N}, psi, i, j)

Compute ⟨ψ|SᶻᵢSᶻⱼ|ψ⟩ in the given basis.
"""
function szsz_correlation(basis::SpinBasis{N}, psi::AbstractVector, i::Int, j::Int) where {N}
    corr = 0.0
    for (idx, state) in enumerate(basis.states)
        _, coeff = apply_SzSz(state, i, j)
        corr += coeff * abs2(psi[idx])
    end
    return real(corr)
end
