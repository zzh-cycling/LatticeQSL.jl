# ═══════════════════════════════════════════════════════════════
# Entanglement observables
# ═══════════════════════════════════════════════════════════════

"""
    reduced_density_matrix(basis::SpinBasis{N}, subsystems::Vector{Int}, state::AbstractVector)

Compute the reduced density matrix of `state` on the subsystem specified by
1-based site indices in `subsystems`.
"""
function reduced_density_matrix(basis::SpinBasis{N}, subsystems::Vector{Int},
                                state::AbstractVector{<:Number}) where {N}
    @assert length(state) == basis.dim "expected state length $(basis.dim), got $(length(state))"

    subsystem = sort!(unique(subsystems))
    @assert !isempty(subsystem) "subsystems must be non-empty"
    @assert all(1 <= site <= N for site in subsystem) "subsystem sites must lie in 1:$N"

    environment = [site for site in 1:N if site ∉ subsystem]
    dim_a = 1 << length(subsystem)
    dim_b = 1 << length(environment)
    psi = zeros(ComplexF64, dim_a, dim_b)

    for (idx, basis_state) in enumerate(basis.states)
        amp = ComplexF64(state[idx])
        abs(amp) < 1e-14 && continue
        ia = _substate_index(basis_state, subsystem)
        ib = _substate_index(basis_state, environment)
        psi[ia + 1, ib + 1] += amp
    end

    return psi * psi'
end

"""
    reduced_density_matrix(basis::SymmetryReducedBasis{N}, subsystems::Vector{Int}, state::AbstractVector)

Expand the reduced-basis state to the parent basis, then trace out the complement.
"""
function reduced_density_matrix(basis::SymmetryReducedBasis{N}, subsystems::Vector{Int},
                                state::AbstractVector{<:Number}) where {N}
    return reduced_density_matrix(basis.parent_basis, subsystems, expand_state(basis, state))
end

"""
    entanglement_entropy(rdm::AbstractMatrix)

Compute the von Neumann entropy `S = -tr(ρ log ρ)`.
"""
function entanglement_entropy(rdm::AbstractMatrix{<:Number})
    @assert ishermitian(rdm) "reduced density matrix must be Hermitian"

    entropy = 0.0
    for value in eigvals(Hermitian(Matrix{ComplexF64}(rdm)))
        p = real(value)
        p <= 1e-12 && continue
        entropy -= p * log(p)
    end
    return entropy
end

"""
    entanglement_entropy(basis, subsystems, state)

Convenience wrapper combining `reduced_density_matrix` and `entanglement_entropy`.
"""
function entanglement_entropy(basis::Union{SpinBasis, SymmetryReducedBasis},
                              subsystems::Vector{Int}, state::AbstractVector{<:Number})
    return entanglement_entropy(reduced_density_matrix(basis, subsystems, state))
end

entanglement_entropy_placeholder() = error("Use entanglement_entropy(...) instead")

function _substate_index(state::BitStr{N, Int}, sites::Vector{Int}) where {N}
    index = 0
    for (offset, site) in enumerate(sites)
        if spin_up(state, site)
            index |= 1 << (offset - 1)
        end
    end
    return index
end
