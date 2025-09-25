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

function _generate_sites(lattice_vectors, lattice_sites, repeats::Vararg{Int,D}; scale = 1.0, verbose::Bool=false) where {D}
    @assert length(lattice_vectors) == D
    @assert D > 0 && length(lattice_sites) > 0
    @assert all(>=(0), repeats)
    @assert all(x -> length(x) == (D), lattice_vectors)
    @assert all(x -> length(x) == (D), lattice_sites)
    T = eltype(first(lattice_vectors))
    locations = NTuple{D,T}[]  # we might want to avoid using `push!` later.
    for ci in CartesianIndices(repeats)
        baseloc = mapreduce(i -> (ci.I[i] - 1) .* lattice_vectors[i], (x, y) -> x .+ y, 1:D)
        for siteloc in lattice_sites
            verbose && @show baseloc, siteloc
            push!(locations, (baseloc .+ siteloc) .* scale)
        end
    end
    return locations
end

"""
    Lattice{D,K,T} <: AbstractLattice{D}
    Lattice(vectors, sites)

The general lattice type for tiling the space, where translation symmetry is assumed. Type parameter `D` is the dimension,
`K` is the number of sites in a unit cell and `T` is the data type for coordinates, e.g. `Float64`. Input arguments are

* `vectors` is a vector/tuple of D-tuple. Its length is D, it specifies the Bravais lattice vectors.
* `sites` is a vector/tuple of D-tuple. Its length is K, it specifies the sites inside a Bravais cell.
* `reciprocal_vectors` is a vector/tuple of D-tuple. Its length is D, it specifies the reciprocal lattice vectors.
* `reciprocal_sites` is a vector/tuple of D-tuple. Its length is K, it specifies the sites inside a reciprocal sites of Bravais cell.
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
        b = ntuple(i -> Tuple(invM[:, i]), Val(D))   # convert to NTuple{D,T}

        τ = 2π
        rsts = ntuple(k -> Tuple(τ .* (invM * collect(sites[k]))), Val(K))
        
        return new{D,K,T}(vectors, sites, b, rsts)
    end
end
_datatype(x::Tuple) = promote_type(_datatype.(x)...)
_datatype(x::Number) = typeof(x)
Lattice(vectors, sites) = Lattice((Tuple.(vectors)...,), (Tuple.(sites)...,))

############ manipulate grid ###############
"""
    MaskedGrid{T}
    MaskedGrid(xs, ys, mask)

Masked square lattice contains 3 fields, the x-coordinates, y-coordinates and a mask, returning true if the site exists, false otherwise. To index the given lattice in the square grid.
e.g. `MaskedGrid([0.0, 1.0, 3.0], [0.0, 2.0,6.0], Bool[1 0 0; 0 1 1; 0 1 0])` specifies the following lattice:

         y₁   y₂        y₃
         ↓    ↓         ↓
    x₁ → ●    ⋅         ●
    x₂ → ⋅    ●         ●

    x₃ → ⋅    ●         ⋅
"""
struct MaskedGrid{T}
    xs::Vector{T}
    ys::Vector{T}
    mask::Matrix{Bool}
end


"""
    make_grid(sites; atol=...)

Create a [`MaskedGrid`](@ref) from the sites. It is required by lattice preparation of Rydberg array.
Because the grid will sort the sites by rows, we need `atol` (default value is 10 time sit data precision)
determines up to what level of round off error, two atoms belong to the same row.

```jldoctest
julia> sites = HoneycombLattice() |> generate_sites(_, 2, 3)
12-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (0.5, 0.2886751345948129)
 (1.0, 0.0)
 (1.5, 0.2886751345948129)
 (0.5, 0.8660254037844386)
 (1.0, 1.1547005383792515)
 (1.5, 0.8660254037844386)
 (2.0, 1.1547005383792515)
 (1.0, 1.7320508075688772)
 (1.5, 2.0207259421636903)
 (2.0, 1.7320508075688772)
 (2.5, 2.0207259421636903)
julia> make_grid(sites)
LatticeQSL.MaskedGrid{Float64}([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], [0.0, 0.2886751345948129, 0.8660254037844386, 1.1547005383792515, 1.7320508075688772, 2.0207259421636903], Bool[1 0 0 0 0 0; 0 1 1 0 0 0; 1 0 0 1 1 0; 0 1 1 0 0 1; 0 0 0 1 1 0; 0 0 0 0 0 1])
``` 

It will generate grid like:
1  0  0  0  0  0
0  1  1  0  0  0
1  0  0  1  1  0
0  1  1  0  0  1
0  0  0  1  1  0
0  0  0  0  0  1
"""
function make_grid(sites::Vector{Tuple{T, T}}; atol = 10 * eps(T)) where {T}
    # sites = padydim(sites)
    xs = sort!(approximate_unique(getindex.(sites, 1), atol))
    ys = sort!(approximate_unique(getindex.(sites, 2), atol))
    ixs = map(s -> findfirst(==(s[1]), xs), sites)
    iys = map(s -> findfirst(==(s[2]), ys), sites)
    m, n = length(xs), length(ys)
    mask = zeros(Bool, m, n)
    for (ix, iy) in zip(ixs, iys)
        mask[ix, iy] = true
    end
    return MaskedGrid(xs, ys, mask)
end
make_grid(sites::Vector{Tuple{T}}; atol = 10 * eps(T)) where {T} = make_grid(map(s -> (s[1], zero(T)), sites); atol = atol)

# return `(uxs, ixs)``, where `uxs` is the unique x-coordinates, `ixs` the mapping from the index in `xs` to the index in `uxs`.
function approximate_unique(xs::AbstractVector{T}, atol) where {T}
    uxs = T[]
    for x in xs
        found = false
        for ux in uxs
            if isapprox(x, ux; atol = atol)
                found = true
                break
            end
        end
        if !found
            push!(uxs, x)
        end
    end
    return uxs
end

function collect_index(mg::MaskedGrid)
    return map(ci -> (mg.xs[ci.I[1]], mg.ys[ci.I[2]]), findall(mg.mask))
end

function grid_index(lt::Vector{NTuple{D, T}}) where {D, T}
    index = Dict{NTuple{D,T}, Int}()
    mask_grid_index = collect_index(make_grid(lt))
    for (i, site) in enumerate(mask_grid_index)
        index[site] = i
    end
    return index
end

"""
    generate_sites(lattice::AbstractLattice{D}, repeats::Vararg{Int,D}; scale=1.0)

Returns an array of tuples (lattice coordinates) by tiling the specified `lattice`.
The tiling repeat the `sites` of the lattice `m` times along the first dimension,
`n` times along the second dimension, and so on. `scale` is a real number that re-scales the lattice constant and atom locations.
"""
function generate_sites(lattice::AbstractLattice{D}, repeats::Vararg{Int,D}; scale = 1.0) where {D}
    return _generate_sites((lattice_vectors(lattice)...,), (lattice_sites(lattice)...,), repeats...; scale = scale)
end

"""
    lattice_vectors(lattice::AbstractLattice)
    reciprocal_lattice_vectors(lattice::AbstractLattice)

Returns Bravais lattice vectors as a D-Tuple of D-Tuple, where D is the space dimension, and its reciprocal lattice vectors.
"""
lattice_vectors(lattice::Lattice) = lattice.vectors
reciprocal_lattice_vectors(lattice::Lattice) = lattice.reciprocal_vectors

"""
    lattice_sites(lattice::AbstractLattice)
    reciprocal_lattice_site(lattice::AbstractLattice)

Returns sites in a Bravais lattice unit cell as a Tuple of D-Tuple, where D is the space dimension, and its reciprocal lattice sites.
"""
lattice_sites(lattice::Lattice) = lattice.sites
reciprocal_lattice_site(lattice::Lattice) = lattice.reciprocal_sites


"""
    struct HoneycombLattice <: AbstractLattice{2}

Type representing 2D Honeycomb Lattice.

Used as an argument of [`generate_sites`](@ref) function to produce tiling in a honeycomb pattern, 
with number of site repetitions being specified by additional arguments of [`generate_sites`](@ref).
Honeycomb is a 2D Lattice, so there must be two integer arguments as additional inputs.

# Example

```jldoctest
julia> generate_sites(HoneycombLattice(), 5, 5)
50-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (0.5, 0.2886751345948129)
 (1.0, 0.0)
 (1.5, 0.2886751345948129)
 (2.0, 0.0)
 (2.5, 0.2886751345948129)
 (3.0, 0.0)
 (3.5, 0.2886751345948129)
 (4.0, 0.0)
 (4.5, 0.2886751345948129)
 ⋮
 (2.5, 3.7527767497325675)
 (3.0, 3.4641016151377544)
 (3.5, 3.7527767497325675)
 (4.0, 3.4641016151377544)
 (4.5, 3.7527767497325675)
 (5.0, 3.4641016151377544)
 (5.5, 3.7527767497325675)
 (6.0, 3.4641016151377544)
 (6.5, 3.7527767497325675)

```
Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::HoneycombLattice)`](@ref) and
[`lattice_sites(::HoneycombLattice)`](@ref).
"""
struct HoneycombLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::HoneycombLattice)
    reciprocal_lattice_vectors(::HoneycombLattice)

Returns the Bravais lattice vectors for a Honeycomb lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice vectors, satisfying 𝐚ᵢ · bⱼ = 2πδᵢⱼ.

The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.5, 0.5√3)
- b₁ = (2π, -2π/√3)
- b₂ = (0.0, 4π/√3)
"""
lattice_vectors(::HoneycombLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))
reciprocal_lattice_vectors(::HoneycombLattice) = ((2π, -2π/sqrt(3)), (0.0, 4π/sqrt(3)))

"""
    lattice_sites(::HoneycombLattice)
    reciprocal_lattice_sites(::HoneycombLattice)

Returns the Bravais Lattice sites for a Honeycomb lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice sites.

The sites are defined as:
- (0.0, 0.0)
- (0.5, 0.5√3)
its reciprocal lattice sites are defined as:
- (0.0, 0.0)
- (1/3, 1/√3)
"""
lattice_sites(::HoneycombLattice) = ((0.0, 0.0), (0.5, 0.5 / sqrt(3)))
reciprocal_lattice_sites(::HoneycombLattice) = ((0.0, 0.0), (1/3, 1/sqrt(3)))

"""
    struct SquareLattice <: AbstractLattice{2}

Type representing 2D Square Lattice.

Used as an argument of [`generate_sites`](@ref) function to produce tiling in a square pattern,
with number of site repetitions being specified by additional arguments of [`generate_sites`](@ref).
Square is a 2D Lattice, so there must be two integer arguments as additional inputs.

# Example

```jldoctest
julia> generate_sites(SquareLattice(), 4, 4)
16-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (1.0, 0.0)
 (2.0, 0.0)
 (3.0, 0.0)
 (0.0, 1.0)
 (1.0, 1.0)
 (2.0, 1.0)
 (3.0, 1.0)
 (0.0, 2.0)
 (1.0, 2.0)
 (2.0, 2.0)
 (3.0, 2.0)
 (0.0, 3.0)
 (1.0, 3.0)
 (2.0, 3.0)
 (3.0, 3.0)

```
Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::SquareLattice)`](@ref) and
[`lattice_sites(::SquareLattice)`](@ref).
"""
struct SquareLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::SquareLattice)
    reciprocal_lattice_vectors(::SquareLattice)

Returns the Bravais lattice vectors for a Square lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice vectors.
    
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.0, 1.0)
- b₁ = (2π, 0.0)
- b₂ = (0.0, 2π)
"""
lattice_vectors(::SquareLattice) = ((1.0, 0.0), (0.0, 1.0))
reciprocal_lattice_vectors(::SquareLattice) = ((2π, 0.0), (0.0, 2π))

"""
    lattice_sites(::SquareLattice)
    reciprocal_lattice_sites(::SquareLattice)

Returns the Bravais Lattice sites for a Square lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice sites.

The sites are defined as:
- (0.0, 0.0)
its reciprocal lattice sites are defined as:
- (0.0, 0.0)
"""
lattice_sites(::SquareLattice) = ((0.0, 0.0),)
reciprocal_lattice_sites(::SquareLattice) = ((0.0, 0.0),)

"""
    struct TriangularLattice <: AbstractLattice{2}

Type representing 2D Square Lattice.

Used as an argument of [`generate_sites`](@ref) function to produce tiling in a triangle pattern, 
with number of site repetitions being specified by additional arguments of [`generate_sites`](@ref).
Triangle is a 2D Lattice, so there must be two integer arguments as additional inputs.

# Example

```jldoctest
julia> generate_sites(TriangularLattice(), 3, 3)
9-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (1.0, 0.0)
 (2.0, 0.0)
 (0.5, 0.8660254037844386)
 (1.5, 0.8660254037844386)
 (2.5, 0.8660254037844386)
 (1.0, 1.7320508075688772)
 (2.0, 1.7320508075688772)
 (3.0, 1.7320508075688772)

```
Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::TriangularLattice)`](@ref) and
[`lattice_sites(::TriangularLattice)`](@ref).
"""
struct TriangularLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::TriangularLattice)

Returns the Bravais lattice vectors for a Triangular lattice as a Tuple of Tuples containing
floats.
    
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.5, 0.5√3)
"""
lattice_vectors(::TriangularLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))

"""
    lattice_sites(::TriangularLattice)

Returns the Bravais Lattice sites for a Triangular lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
"""
lattice_sites(::TriangularLattice) = ((0.0, 0.0),)

"""
    struct ChainLattice <: AbstractLattice{1}

Type representing 1D Chain Lattice.

Used as an argument of [`generate_sites`](@ref) function to produce tiling in a chain pattern, 
with number of site repetitions being specified by the additional argument of [`generate_sites`](@ref).
Chain is a 1D Lattice, so there must be one integer argument as additional inputs.

# Example

```jldoctest
julia> generate_sites(ChainLattice(), 5)
5-element Vector{Tuple{Float64}}:
 (0.0,)
 (1.0,)
 (2.0,)
 (3.0,)
 (4.0,)
```
Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::ChainLattice)`](@ref) and
[`lattice_sites(::ChainLattice)`](@ref).
"""
struct ChainLattice <: AbstractLattice{1} end

"""
    lattice_vectors(::ChainLattice)
    reciprocal_lattice_vectors(::ChainLattice)

Returns the Bravais lattice vectors for a Chain lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice vectors.
    
The vectors are defined as:
- 𝐚₁ = (1.0,)
- b₁ = (2π,)
"""
lattice_vectors(::ChainLattice) = ((1.0,),)
reciprocal_lattice_vectors(::ChainLattice) = ((2π,),)

"""
    lattice_sites(::ChainLattice)
    reciprocal_lattice_sites(::ChainLattice)

Returns the Bravais Lattice sites for a Chain lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0,)
its reciprocal lattice sites are defined as:
- (0.0,)
"""
lattice_sites(::ChainLattice) = ((0.0,),)
reciprocal_lattice_sites(::ChainLattice) = ((0.0,),)

"""
    struct LiebLattice <: AbstractLattice{2}

Type representing 2D Lieb Lattice.

Used as an argument of [`generate_sites`](@ref) function to produce tiling in a Lieb (square-depleted) pattern, 
with number of site repetitions being specified by additional arguments of [`generate_sites`](@ref).
Lieb is a 2D Lattice, so there must be two integer arguments as additional inputs.

# Example

```jldoctest
julia> generate_sites(LiebLattice(), 3, 3)
27-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (0.5, 0.0)
 (0.0, 0.5)
 (1.0, 0.0)
 (1.5, 0.0)
 (1.0, 0.5)
 (2.0, 0.0)
 (2.5, 0.0)
 (2.0, 0.5)
 (0.0, 1.0)
 ⋮
 (0.0, 2.0)
 (0.5, 2.0)
 (0.0, 2.5)
 (1.0, 2.0)
 (1.5, 2.0)
 (1.0, 2.5)
 (2.0, 2.0)
 (2.5, 2.0)
 (2.0, 2.5)
```
Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::LiebLattice)`](@ref) and
[`lattice_sites(::LiebLattice)`](@ref).
"""
struct LiebLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::LiebLattice)
    reciprocal_lattice_vectors(::LiebLattice)

Returns the Bravais lattice vectors for a Lieb lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice vectors.
        
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.0, 1.0)
- b₁ = (2π, 0.0)
- b₂ = (0.0, 2π)
"""
lattice_vectors(::LiebLattice) = ((1.0, 0.0), (0.0, 1.0))
reciprocal_lattice_vectors(::LiebLattice) = ((2π, 0.0), (0.0, 2π))

"""
    lattice_sites(::LiebLattice)
    reciprocal_lattice_sites(::LiebLattice)

Returns the Bravais Lattice sites for a Lieb lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
- (0.5, 0.0)
- (0.0, 0.5)
its reciprocal lattice sites are defined as:
- (0.0, 0.0)

"""
lattice_sites(::LiebLattice) = ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5))

"""
    struct KagomeLattice <: AbstractLattice{2}

Type representing 2D Kagome Lattice.

Used as an argument of [`generate_sites`](@ref) function to produce tiling in a Kagome pattern, 
with number of site repetitions being specified by additional arguments of [`generate_sites`](@ref).
Kagome is a 2D Lattice, so there must be two integer arguments as additional inputs.

# Example

```jldoctest
julia> generate_sites(KagomeLattice(), 3, 3)
27-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (0.25, 0.4330127018922193)
 (0.75, 0.4330127018922193)
 (1.0, 0.0)
 (1.25, 0.4330127018922193)
 (1.75, 0.4330127018922193)
 (2.0, 0.0)
 (2.25, 0.4330127018922193)
 (2.75, 0.4330127018922193)
 (0.5, 0.8660254037844386)
 ⋮
 (1.0, 1.7320508075688772)
 (1.25, 2.1650635094610964)
 (1.75, 2.1650635094610964)
 (2.0, 1.7320508075688772)
 (2.25, 2.1650635094610964)
 (2.75, 2.1650635094610964)
 (3.0, 1.7320508075688772)
 (3.25, 2.1650635094610964)
 (3.75, 2.1650635094610964)
```
Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::KagomeLattice)`](@ref) and
[`lattice_sites(::KagomeLattice)`](@ref).
"""
struct KagomeLattice <: AbstractLattice{2} end

"""
    lattice_vectors(::KagomeLattice)
    reciprocal_lattice_vectors(::KagomeLattice)

Returns the Bravais lattice vectors for a Kagome lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice vectors.
        
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.5, 0.5√3)
- b₁ = (2π, -2π/√3)
- b₂ = (0.0, 4π/√3)
"""
lattice_vectors(::KagomeLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))
reciprocal_lattice_vectors(::KagomeLattice) = ((2π, -2π/√3), (0.0, 4π/√3))

"""
    lattice_sites(::KagomeLattice)
    reciprocal_lattice_sites(::KagomeLattice)

Returns the Bravais Lattice sites for a Kagome lattice as a Tuple of Tuples containing
floats.

The sites are defined as:
- (0.0, 0.0)
- (0.25, 0.25√3)
- (0.75, 0.25√3)
its reciprocal lattice sites are defined as:
- (0.0, 0.0)
- (1/3, 1/√3)
- (2/3, 1/√3)
"""
lattice_sites(::KagomeLattice) = ((0.0, 0.0), (0.25, 0.25 * sqrt(3)), (0.75, 0.25 * sqrt(3)))
reciprocal_lattice_sites(::KagomeLattice) = ((0.0, 0.0), (1/3, 1/sqrt(3)), (2/3, 1/sqrt(3)))

"""
    struct RectangularLattice <: AbstractLattice{2}

Type representing 2D Rectangular Lattice.

Used as an argument of [`generate_sites`](@ref) function to produce tiling in a Rectangular pattern, 
with number of site repetitions being specified by additional arguments of [`generate_sites`](@ref).
Rectangular is a 2D Lattice, so there must be two integer arguments as additional inputs.
This type also enables the user to modify the length of one of the Bravais lattice vectors, 
by passing a single integer on construction as the aspect ratio.

# Example

```jldoctest
julia> generate_sites(RectangularLattice(2.0), 2, 2)
4-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (1.0, 0.0)
 (0.0, 2.0)
 (1.0, 2.0)
```
Overriden functions to return lattice vectors and sites exists as 
[`lattice_vectors(::RectangularLattice)`](@ref) and
[`lattice_sites(::RectangularLattice)`](@ref).
"""
struct RectangularLattice <: AbstractLattice{2}
    aspect_ratio::Float64
end

"""
    lattice_vectors(r::RectangularLattice)
    reciprocal_lattice_vectors(r::RectangularLattice)

Returns the Bravais lattice vectors for a Rectangular lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice vectors.
        
The vectors are defined as:
- 𝐚₁ = (1.0, 0.0)
- 𝐚₂ = (0.0, `r.aspect_ratio`), where `aspect_ratio` is a `Float64`.
- b₁ = (2π, 0.0)
- b₂ = (0.0, 2π/`r.aspect_ratio`)
"""
lattice_vectors(r::RectangularLattice) = ((1.0, 0.0), (0.0, r.aspect_ratio))
reciprocal_lattice_vectors(r::RectangularLattice) = ((2π, 0.0), (0.0, 2π / r.aspect_ratio))

"""
    lattice_sites(::RectangularLattice)
    reciprocal_lattice_sites(::RectangularLattice)

Returns the Bravais Lattice sites for a Rectangular lattice as a Tuple of Tuples containing
floats, and its reciprocal lattice sites.

The sites are defined as:
- (0.0, 0.0)
its reciprocal lattice sites are defined as:
- (0.0, 0.0)
"""
lattice_sites(::RectangularLattice) = ((0.0, 0.0),)
reciprocal_lattice_sites(::RectangularLattice) = ((0.0, 0.0),)


############ manipulate sites ###############
"""
    offset_axes(sites, offset0, offset1, ...)

Offset the `sites` by distance specified by `offset0`, `offset1`, ...

```jldoctest
julia> sites = [(1.0, 2.0), (10.0, 3.0), (1.0, 12.0), (3.0, 5.0)]
4-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
 (10.0, 3.0)
 (1.0, 12.0)
 (3.0, 5.0)

julia> LatticeQSL.offset_axes(sites, 1.0, 3.0)
4-element Vector{Tuple{Float64, Float64}}:
 (2.0, 5.0)
 (11.0, 6.0)
 (2.0, 15.0)
 (4.0, 8.0)
```
"""
function offset_axes(sites, offset0::T, offsets::Vararg{T,D}) where {D,T}
    @assert all(x -> length(x) == D + 1, sites) "expected $(D + 1)-tuple sites, got $(length.(sites))"
    return map(x -> (x[1] + offset0, (x[2:end] .+ offsets)...), sites)
end

"""
    rescale_axes(sites, scale::Real)

Rescale the `sites` by a constant `scale`.

```jldoctest
julia> sites = [(1.0, 2.0), (10.0, 3.0), (1.0, 12.0), (3.0, 5.0)]
4-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
 (10.0, 3.0)
 (1.0, 12.0)
 (3.0, 5.0)

julia> LatticeQSL.rescale_axes(sites, 2.0)
4-element Vector{Tuple{Float64, Float64}}:
 (2.0, 4.0)
 (20.0, 6.0)
 (2.0, 24.0)
 (6.0, 10.0)
```
"""
function rescale_axes(sites, scale::Real)
    return map(x -> x .* scale, sites)
end

"""
    random_dropout(sites, ratio::Real)

Randomly drop out `ratio * number of sites` atoms from `sites`, where `ratio` ∈ [0, 1].
"""
function random_dropout(sites::AbstractVector, ratio::Real)
    0 ≤ ratio ≤ 1 || throw(ArgumentError("ratio must be ∈ [0, 1], got $ratio"))
    k = floor(Int, length(sites) * (1 - ratio))
    k < 1 && return empty(sites)          # delete all
    k ≥ length(sites) && return copy(sites)  # keep all
    ind = randperm(length(sites))[1:k]  # first k unique random indices
    sites[sort!(ind)]                   # keep original order
end

"""
    clip_axes(sites, bound0, bound1, ...)

Remove sites out of `bounds`, where `bounds` is specified by `bound0`, `bound1`, ...

```jldoctest
julia> sites = [(1.0, 2.0), (10.0, 3.0), (1.0, 12.0), (3.0, 5.0)]
4-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
 (10.0, 3.0)
 (1.0, 12.0)
 (3.0, 5.0)

julia> LatticeQSL.clip_axes(sites, (-5.0, 5.0), (-5.0, 5.0))
2-element Vector{Tuple{Float64, Float64}}:
 (1.0, 2.0)
 (3.0, 5.0)
```
"""
function clip_axes(sites, bound0::Tuple{T,T}, bounds::Vararg{Tuple{T,T},D}) where {D,T}
    @assert all(x -> length(x) == D + 1, sites) "expected $(D + 1)-tuple sites, got $(length.(sites))"
    return filter(x -> bound0[1] <= x[1] <= bound0[2] && all(i -> bounds[i][1] <= x[i+1] <= bounds[i][2], 1:D), sites)
end