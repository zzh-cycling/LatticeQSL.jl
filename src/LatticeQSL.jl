module LatticeQSL

using LinearAlgebra, SparseArrays, Random, Arpack, BitBasis

export Lattice, lattice_vectors, generate_sites, AbstractLattice,
    HoneycombLattice,
    SquareLattice,
    TriangularLattice,
    ChainLattice,
    LiebLattice,
    KagomeLattice,
    GeneralLattice,
    RectangularLattice,
    # interfaces
    generate_sites,
    offset_axes,
    random_dropout,
    rescale_axes,
    clip_axes,
    lattice_sites,
    lattice_vectors,
    make_grid
    
export kitaev_hamiltonian_sparse, wilson12, loop_path, loop_op, honeycomb_strings
export honeycomb_basis, kitaev_honeycomb_ham
export kagome_strings, kagome_basis, Heisenberg_hamiltonian_sparse, loop_map

include("Lattice.jl")
include("EDdemo.jl")
include("Basis.jl")
end
