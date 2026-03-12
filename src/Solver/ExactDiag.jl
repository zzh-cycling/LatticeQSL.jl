# ═══════════════════════════════════════════════════════════════
# Exact diagonalization utilities
# ═══════════════════════════════════════════════════════════════

using LinearAlgebra

"""
    exact_diagonalization(H::AbstractMatrix)

Compute all eigenvalues and eigenvectors of the Hamiltonian `H`.
For small systems only.

Returns `(eigenvalues, eigenvectors)` sorted in ascending energy.
"""
function exact_diagonalization(H::AbstractMatrix)
    F = eigen(Matrix(H))
    return F.values, F.vectors
end

"""
    ground_state_exact(H::AbstractMatrix; per_site::Bool=false, nsites::Int=1)

Compute the ground-state energy and vector using full diagonalization.
"""
function ground_state_exact(H::AbstractMatrix; per_site::Bool = false, nsites::Int = 1)
    vals, vecs = exact_diagonalization(H)
    e0 = vals[1]
    per_site && (e0 /= nsites)
    return e0, vecs[:, 1]
end
