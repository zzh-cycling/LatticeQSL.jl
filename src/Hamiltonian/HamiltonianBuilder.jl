# ═══════════════════════════════════════════════════════════════
# HamiltonianBuilder: build sparse Hamiltonians from lattices, bases, and models
#
# No Kronecker products are used.
# We iterate over basis states, apply local terms, and fill a SparseMatrixCSC.
# ═══════════════════════════════════════════════════════════════

using SparseArrays

"""
    build_hamiltonian(fl::FiniteLattice, basis::SpinBasis, model::AbstractModel)

Construct the sparse Hamiltonian matrix in the given basis.
Returns a `SparseMatrixCSC` of size `dim × dim`.
"""
function build_hamiltonian(fl::FiniteLattice, basis::SpinBasis{N}, model::AbstractModel) where {N}
    row_inds = Int[]
    col_inds = Int[]
    values = ComplexF64[]
    
    bonds = model_bonds(model, fl)
    
    for col in 1:basis.dim
        state = basis[col]
        
        for (i, j, coupling, bond_type) in bonds
            results = apply_model_bond(model, state, i, j, bond_type)
            for (new_state, coeff) in results
                row = state_index(basis, new_state)
                row == 0 && continue
                push!(row_inds, row)
                push!(col_inds, col)
                push!(values, coupling * ComplexF64(coeff))
            end
        end
        
        # Add transverse field for Ising model: h Σᵢ Sˣᵢ
        if model isa IsingModel && model.h != 0
            for site in 1:fl.nsites
                new_state, coeff = apply_Sx(state, site)
                row = state_index(basis, new_state)
                row == 0 && continue  # may leave Sz sector
                push!(row_inds, row)
                push!(col_inds, col)
                push!(values, model.h * ComplexF64(coeff))
            end
        end
    end
    
    H = sparse(row_inds, col_inds, values, basis.dim, basis.dim)
    return dropzeros!(H)
end

"""
    build_hamiltonian(fl::FiniteLattice, basis::SymmetryReducedBasis, model::AbstractModel)

Construct the sparse Hamiltonian matrix in a symmetry-reduced basis.
"""
function build_hamiltonian(fl::FiniteLattice, basis::SymmetryReducedBasis{N}, model::AbstractModel) where {N}
    row_inds = Int[]
    col_inds = Int[]
    values = ComplexF64[]
    
    bonds = model_bonds(model, fl)
    
    for col in 1:basis.dim
        rep_state = basis[col]
        
        for (i, j, coupling, bond_type) in bonds
            results = apply_model_bond(model, rep_state, i, j, bond_type)
            for (new_state, coeff) in results
                row, factor = matrix_element_factor(basis, col, new_state)
                row == 0 && continue
                push!(row_inds, row)
                push!(col_inds, col)
                push!(values, coupling * ComplexF64(coeff) * factor)
            end
        end
        
        if model isa IsingModel && model.h != 0
            for site in 1:fl.nsites
                new_state, coeff = apply_Sx(rep_state, site)
                row, factor = matrix_element_factor(basis, col, new_state)
                row == 0 && continue
                push!(row_inds, row)
                push!(col_inds, col)
                push!(values, model.h * ComplexF64(coeff) * factor)
            end
        end
    end
    
    H = sparse(row_inds, col_inds, values, basis.dim, basis.dim)
    return dropzeros!(H)
end

"""
    hamiltonian_dimension(basis) -> Int

Return the Hilbert space dimension for either a full basis or reduced basis.
"""
hamiltonian_dimension(basis::SpinBasis) = basis.dim
hamiltonian_dimension(basis::SymmetryReducedBasis) = basis.dim

"""
    basis_summary(basis)

Return a named tuple summarizing the basis dimensions.
"""
function basis_summary(basis::SpinBasis)
    return (nsites = basis.nsites, dim = basis.dim, reduced = false)
end

function basis_summary(basis::SymmetryReducedBasis)
    return (
        nsites = basis.nsites,
        dim = basis.dim,
        parent_dim = basis.parent_basis.dim,
        reduction_ratio = basis.dim / basis.parent_basis.dim,
        reduced = true,
    )
end
