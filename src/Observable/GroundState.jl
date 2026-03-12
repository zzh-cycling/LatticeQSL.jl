# ═══════════════════════════════════════════════════════════════
# Ground-state observables
# ═══════════════════════════════════════════════════════════════

using LinearAlgebra

"""
    expectation(psi, O)

Compute ⟨ψ|O|ψ⟩.
"""
expectation(psi::AbstractVector, O::AbstractMatrix) = dot(psi, O * psi)

"""
    energy_per_site(e0, nsites)

Convert total ground-state energy to energy per site.
"""
energy_per_site(e0::Real, nsites::Int) = e0 / nsites
