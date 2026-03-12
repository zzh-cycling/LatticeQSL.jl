module LatticeQSL

using LinearAlgebra, SparseArrays, Random, Arpack, BitBasis

include("Lattice/AbstractLattice.jl")
include("Lattice/PredefinedLattices.jl")
include("Lattice/LatticeUtils.jl")
include("Lattice/FiniteLattice.jl")

include("Basis/SpinBasis.jl")
include("Basis/Symmetry.jl")
include("Basis/SpinSymmetry.jl")
include("Basis/TranslationSymmetry.jl")
include("Basis/PointGroupSymmetry.jl")
include("Basis/SymmetryReducedBasis.jl")

include("Hamiltonian/SpinOperators.jl")
include("Hamiltonian/AbstractModel.jl")
include("Hamiltonian/Models.jl")
include("Hamiltonian/HamiltonianBuilder.jl")

include("Solver/ExactDiag.jl")
include("Solver/Lanczos.jl")

include("Observable/GroundState.jl")
include("Observable/Correlation.jl")
include("Observable/StructureFactor.jl")
include("Observable/Entanglement.jl")
include("Observable/WilsonLoop.jl")

# compatibility layer from the current codebase
include("EDdemo.jl")
include("Basis.jl")

export AbstractLattice, dimension, nsites_per_cell,
    Lattice, GeneralLattice,
    HoneycombLattice, SquareLattice, TriangularLattice, ChainLattice,
    KagomeLattice, LiebLattice, RectangularLattice,
    lattice_vectors, reciprocal_lattice_vector, lattice_sites, reciprocal_lattice_sites,
    generate_sites, offset_axes, rescale_axes, random_dropout, clip_axes,
    AtomList, MaskedGrid, make_grid, collect_atoms,
    Bond, FiniteLattice, site_index, neighbors, bonds_by_label, all_bond_labels,
    translation_permutations, translation_group

export SpinBasis, state_index, spin_up, spin_value, flip_spin, total_sz, num_up, permute_bits,
    AbstractSymmetry, AbstractSpatialSymmetry, AbstractInternalSymmetry,
    SzConservation, target_sz, SpinInversion, invert_spins,
    TranslationSymmetry, momentum_values,
    PointGroupSymmetry,
    SymmetryReducedBasis, representative_index, matrix_element_factor, expand_state

export apply_Sz, apply_Sp, apply_Sm, apply_Sx, apply_Sy,
    apply_SzSz, apply_SpSm, apply_SmSp, apply_heisenberg, apply_kitaev_bond,
    AbstractModel, HeisenbergModel, KitaevModel, IsingModel, XXZModel, GeneralSpinModel,
    model_bonds, apply_model_bond,
    build_hamiltonian, hamiltonian_dimension, basis_summary

export exact_diagonalization, ground_state_exact, ground_state, low_energy_spectrum,
    expectation, energy_per_site, szsz_correlation,
    structure_factor_placeholder, reduced_density_matrix, entanglement_entropy,
    entanglement_entropy_placeholder, wilson_loop_placeholder

export kitaev_hamiltonian_sparse, wilson12, loop_path, loop_op, honeycomb_strings,
    honeycomb_basis, kitaev_honeycomb_ham,
    kagome_strings, kagome_basis, Heisenberg_hamiltonian_sparse, loop_map

end
