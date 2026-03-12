# ═══════════════════════════════════════════════════════════════
# Lattice utility functions: generate_sites, MaskedGrid, AtomList, etc.
# ═══════════════════════════════════════════════════════════════

# ─────────── Site generation ───────────

function _generate_sites(lattice_vectors, lattice_sites, repeats::Vararg{Int,D}; scale = 1.0) where {D}
    @assert length(lattice_vectors) == D
    @assert D > 0 && length(lattice_sites) > 0
    @assert all(>=(0), repeats)
    @assert all(x -> length(x) == D, lattice_vectors)
    @assert all(x -> length(x) == D, lattice_sites)
    T = eltype(first(lattice_vectors))
    locations = NTuple{D,T}[]
    for ci in CartesianIndices(repeats)
        baseloc = mapreduce(i -> (ci.I[i] - 1) .* lattice_vectors[i], (x, y) -> x .+ y, 1:D)
        for siteloc in lattice_sites
            push!(locations, (baseloc .+ siteloc) .* scale)
        end
    end
    return locations
end

"""
    generate_sites(lattice::AbstractLattice{D}, repeats::Vararg{Int,D}; scale=1.0)

Returns an array of tuples (lattice coordinates) by tiling the specified `lattice`.
"""
function generate_sites(lattice::AbstractLattice{D}, repeats::Vararg{Int,D}; scale = 1.0) where {D}
    return _generate_sites((lattice_vectors(lattice)...,), (lattice_sites(lattice)...,), repeats...; scale = scale)
end

# ─────────── AtomList ───────────

"""
    AtomList{D,T}

A wrapper around `Vector{NTuple{D,T}}` for convenient atom position manipulation.
"""
struct AtomList{D,T} <: AbstractVector{NTuple{D,T}}
    positions::Vector{NTuple{D,T}}
end
AtomList(v::Vector{<:Tuple}) = AtomList(collect(v))

Base.size(al::AtomList) = size(al.positions)
Base.getindex(al::AtomList, i::Int) = al.positions[i]
Base.getindex(al::AtomList, i::AbstractVector) = AtomList(al.positions[i])
Base.setindex!(al::AtomList, v, i::Int) = (al.positions[i] = v)
Base.length(al::AtomList) = length(al.positions)
Base.iterate(al::AtomList, args...) = iterate(al.positions, args...)
Base.:(==)(a::AtomList, b::AtomList) = a.positions == b.positions
Base.:(==)(a::AtomList, b::Vector) = a.positions == b
Base.:(==)(a::Vector, b::AtomList) = a == b.positions

function Base.deleteat!(al::AtomList, indices::Int...)
    sorted = sort(collect(indices); rev = true)
    for i in sorted
        deleteat!(al.positions, i)
    end
    return al
end
Base.deleteat!(al::AtomList, i::Int) = (deleteat!(al.positions, i); al)

# ─────────── Axes manipulation ───────────

"""
    offset_axes(sites, offset0, offset1, ...)

Offset the `sites` by distance specified by `offset0`, `offset1`, ...
"""
function offset_axes(sites, offset0::T, offsets::Vararg{T,D}) where {D,T}
    @assert all(x -> length(x) == D + 1, sites) "expected $(D + 1)-tuple sites, got $(length.(sites))"
    return map(x -> (x[1] + offset0, (x[2:end] .+ offsets)...), sites)
end

# Curried version for piping
offset_axes(offset0::T, offsets::Vararg{T,D}) where {D,T} = sites -> offset_axes(sites, offset0, offsets...)

"""
    rescale_axes(sites, scale::Real)

Rescale the `sites` by a constant `scale`.
"""
function rescale_axes(sites, scale::Real)
    return map(x -> x .* scale, sites)
end

rescale_axes(scale::Real) = sites -> rescale_axes(sites, scale)

"""
    random_dropout(sites, ratio::Real)

Randomly drop out `ratio * number of sites` atoms from `sites`, where `ratio` ∈ [0, 1].
"""
function random_dropout(sites::AbstractVector, ratio::Real)
    0 ≤ ratio ≤ 1 || throw(ArgumentError("ratio must be ∈ [0, 1], got $ratio"))
    k = floor(Int, length(sites) * (1 - ratio))
    k < 1 && return empty(sites)
    k ≥ length(sites) && return copy(sites)
    ind = randperm(length(sites))[1:k]
    sites[sort!(ind)]
end

random_dropout(ratio::Real) = sites -> random_dropout(sites, ratio)

"""
    clip_axes(sites, bound0, bound1, ...)

Remove sites out of `bounds`, where `bounds` is specified by `bound0`, `bound1`, ...
"""
function clip_axes(sites, bound0::Tuple{T,T}, bounds::Vararg{Tuple{T,T},D}) where {D,T}
    @assert all(x -> length(x) == D + 1, sites) "expected $(D + 1)-tuple sites, got $(length.(sites))"
    return filter(x -> bound0[1] <= x[1] <= bound0[2] && all(i -> bounds[i][1] <= x[i+1] <= bounds[i][2], 1:D), sites)
end

# ─────────── MaskedGrid ───────────

"""
    MaskedGrid{T}

Masked square lattice grid for representing lattice sites on a regular grid.
"""
struct MaskedGrid{T}
    xs::Vector{T}
    ys::Vector{T}
    mask::Matrix{Bool}
end

"""
    collect_atoms(mg::MaskedGrid)

Collect all atom positions from a MaskedGrid.
"""
function collect_atoms(mg::MaskedGrid{T}) where T
    positions = Tuple{T,T}[]
    for ix in 1:length(mg.xs), iy in 1:length(mg.ys)
        if mg.mask[ix, iy]
            push!(positions, (mg.xs[ix], mg.ys[iy]))
        end
    end
    return positions
end

"""
    make_grid(sites; atol=...)

Create a [`MaskedGrid`](@ref) from the sites.
"""
function make_grid(sites::AbstractVector{<:Tuple}; atol = 10 * eps(Float64))
    xs = sort!(approximate_unique(getindex.(sites, 1), atol))
    ys = sort!(approximate_unique(getindex.(sites, 2), atol))
    ixs = map(s -> findfirst(x -> isapprox(x, s[1]; atol), xs), sites)
    iys = map(s -> findfirst(y -> isapprox(y, s[2]; atol), ys), sites)
    m, n = length(xs), length(ys)
    mask = zeros(Bool, m, n)
    for (ix, iy) in zip(ixs, iys)
        mask[ix, iy] = true
    end
    return MaskedGrid(xs, ys, mask)
end

# Handle 1D sites
function make_grid(sites::AbstractVector{Tuple{T}}) where T
    xs = sort!(approximate_unique(getindex.(sites, 1), 10 * eps(T)))
    ys = T[zero(T)]
    mask = ones(Bool, length(xs), 1)
    return MaskedGrid(xs, ys, mask)
end

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
