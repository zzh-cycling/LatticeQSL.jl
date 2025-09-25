module LatticeQSL

using LinearAlgebra, SparseArrays, Random, Arpack, BitBasis

export Lattice, AbstractLattice, RectangularLattice, ChainLattice, SquareLattice, 
    HoneycombLattice, TriangularLattice, LiebLattice, KagomeLattice, 
    # interfaces
    generate_sites,
    offset_axes, clip_axes, random_dropout,
    lattice_sites, lattice_vectors, reciprocal_lattice_sites, reciprocal_lattice_vectors,
    make_grid, collect_index, grid_index
    
export kitaev_hamiltonian_sparse, wilson12, loop_path, loop_op, honeycomb_strings
export honeycomb_basis, kitaev_honeycomb_ham
export kagome_strings, kagome_basis, Heisenberg_hamiltonian_sparse, loop_map

include("Lattice.jl")
include("EDdemo.jl")
include("Basis.jl")
end
