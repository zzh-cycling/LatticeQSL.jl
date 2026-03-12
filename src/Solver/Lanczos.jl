# ═══════════════════════════════════════════════════════════════
# Lanczos / Arnoldi utilities using Arpack
# ═══════════════════════════════════════════════════════════════

using Arpack

"""
    ground_state(H::AbstractMatrix; nev=1, ncv=20, per_site=false, nsites=1)

Find the ground state using sparse diagonalization.
Returns `(e0, psi0)`.
"""
function ground_state(H::AbstractMatrix; nev::Int = 1, ncv::Int = 20,
                      per_site::Bool = false, nsites::Int = 1)
    vals, vecs = eigs(H; nev = nev, which = :SR, ncv = ncv)
    order = sortperm(real.(vals))
    e0 = real(vals[order[1]])
    per_site && (e0 /= nsites)
    return e0, vecs[:, order[1]]
end

"""
    low_energy_spectrum(H::AbstractMatrix; nev=6, ncv=30)

Return the lowest `nev` eigenvalues.
"""
function low_energy_spectrum(H::AbstractMatrix; nev::Int = 6, ncv::Int = 30)
    vals, _ = eigs(H; nev = nev, which = :SR, ncv = ncv)
    return sort(real.(vals))
end
