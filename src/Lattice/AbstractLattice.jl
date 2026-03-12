"""
    AbstractLattice{D}

Supertype for all `D` dimensional lattices.

# Implementation

[`lattice_vectors`](@ref) and [`lattice_sites`](@ref) functions must be defined
which should both return an indexable iterable containing the Bravais lattice vectors and 
lattice sites respectively. (e.g.: [`Lattice`](@ref) returns a tuple of tuples containing the 
Bravais lattice vectors and lattice sites).
"""
abstract type AbstractLattice{D} end

"""
    dimension(::AbstractLattice{D})

Returns the space dimension of target lattice.
e.g. [`ChainLattice`](@ref) is a 1D lattice, hence returns 1.
"""
dimension(::AbstractLattice{D}) where {D} = D

"""
    nsites_per_cell(lattice::AbstractLattice)

Returns the number of sites in the unit cell of the lattice.
"""
nsites_per_cell(lattice::AbstractLattice) = length(lattice_sites(lattice))

"""
    lattice_vectors(lattice::AbstractLattice)

Returns Bravais lattice vectors as a D-Tuple of D-Tuple, where D is the space dimension.
"""
function lattice_vectors end

"""
    reciprocal_lattice_vector(lattice::AbstractLattice)

Returns reciprocal lattice vectors satisfying 𝐚ᵢ · bⱼ = 2πδᵢⱼ.
"""
function reciprocal_lattice_vector end

"""
    lattice_sites(lattice::AbstractLattice)

Returns sites in a Bravais lattice unit cell as a Tuple of D-Tuple, where D is the space dimension.
"""
function lattice_sites end

"""
    reciprocal_lattice_sites(lattice::AbstractLattice)

Returns reciprocal lattice sites.
"""
function reciprocal_lattice_sites end

# ─────────── Concrete Lattice type ───────────

"""
    Lattice{D,K,T} <: AbstractLattice{D}
    Lattice(vectors, sites)

The general lattice type for tiling the space, where translation symmetry is assumed. Type parameter `D` is the dimension,
`K` is the number of sites in a unit cell and `T` is the data type for coordinates, e.g. `Float64`. Input arguments are

* `vectors` is a vector/tuple of D-tuple. Its length is D, it specifies the Bravais lattice vectors.
* `sites` is a vector/tuple of D-tuple. Its length is K, it specifies the sites inside a Bravais cell.
* `reciprocal_vectors` is a vector/tuple of D-tuple. Its length is D, it specifies the reciprocal lattice vectors.
* `reciprocal_sites` is a vector/tuple of D-tuple. Its length is K, it specifies the reciprocal sites of Bravais cell.
"""
struct Lattice{D,K,T} <: AbstractLattice{D}
    vectors::NTuple{D,NTuple{D,T}}
    sites::NTuple{K,NTuple{D,T}}
    reciprocal_vectors::NTuple{D,NTuple{D,T}}
    reciprocal_sites::NTuple{K,NTuple{D,T}}
    function Lattice(vectors::NTuple{D}, sites::NTuple{K}) where {D,K}
        if D == 0 || K == 0
            error("Lattice requires at least one vector and one site")
        end
        T = promote_type(_datatype(vectors), _datatype(sites))

        M = reduce(hcat, [collect(v) for v in vectors])   # basis matrix
        invM = inv(M)'
        b = ntuple(i -> Tuple(invM[:, i]), Val(D))        # reciprocal vectors

        τ = 2π
        rsts = ntuple(k -> Tuple(τ .* (invM * collect(sites[k]))), Val(K))

        return new{D,K,T}(vectors, sites, b, rsts)
    end
end
_datatype(x::Tuple) = promote_type(_datatype.(x)...)
_datatype(x::Number) = typeof(x)
Lattice(vectors, sites) = Lattice((Tuple.(vectors)...,), (Tuple.(sites)...,))

lattice_vectors(lattice::Lattice) = lattice.vectors
reciprocal_lattice_vector(lattice::Lattice) = lattice.reciprocal_vectors
lattice_sites(lattice::Lattice) = lattice.sites
reciprocal_lattice_sites(lattice::Lattice) = lattice.reciprocal_sites
