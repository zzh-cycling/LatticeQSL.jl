# ═══════════════════════════════════════════════════════════════
# FiniteLattice: a finite cluster with boundary conditions,
# neighbor tables, and symmetry permutations for ED
# ═══════════════════════════════════════════════════════════════

"""
    Bond

A bond connecting site `i` to site `j` with a label (e.g. `:x`, `:y`, `:z` for Kitaev,
or `:nn`, `:nnn` for nearest/next-nearest neighbor).
"""
struct Bond
    i::Int
    j::Int
    label::Symbol
end

"""
    FiniteLattice{D, L <: AbstractLattice{D}}

A finite cluster built from a unit cell lattice type `L`, tiled `size` times along each
Bravais direction, with periodic or open boundary conditions.

# Fields
- `unitcell::L` — the infinite lattice unit cell definition
- `size::NTuple{D, Int}` — number of unit cells along each direction, e.g. `(Lx, Ly)`
- `pbc::NTuple{D, Bool}` — periodic boundary conditions per direction
- `nsites::Int` — total number of sites
- `site_positions::Vector{NTuple{D, Float64}}` — Cartesian position of each site
- `bonds::Vector{Bond}` — all bonds (nearest neighbor by default)
- `cell_indices::Array{Int, D}` — maps unit cell grid index `(i₁,...,iD)` → first site index in that cell (1-based)
- `sublattice::Vector{Int}` — sublattice index for each site (1-based)

The site ordering convention is: site index `s` for unit cell at grid position `(c₁,...,cD)`
with sublattice `α` is `s = (linear_index(c) - 1) * K + α`, where `K = nsites_per_cell`.
"""
struct FiniteLattice{D, L <: AbstractLattice{D}}
    unitcell::L
    size::NTuple{D, Int}
    pbc::NTuple{D, Bool}
    nsites::Int
    site_positions::Vector{NTuple{D, Float64}}
    bonds::Vector{Bond}
    cell_indices::Array{Int, D}  # cell grid → first site index of that cell
    sublattice::Vector{Int}
end

"""
    FiniteLattice(unitcell::AbstractLattice{D}, size::NTuple{D,Int}; pbc=ntuple(_->true, D))

Construct a finite lattice cluster from a unit cell definition.
"""
function FiniteLattice(unitcell::AbstractLattice{D}, sz::NTuple{D,Int};
                       pbc::NTuple{D,Bool} = ntuple(_ -> true, Val(D))) where {D}
    vecs = lattice_vectors(unitcell)
    sites = lattice_sites(unitcell)
    K = length(sites)
    ncells = prod(sz)
    nsites = ncells * K

    # Build site positions and sublattice labels
    site_positions = Vector{NTuple{D, Float64}}(undef, nsites)
    sublat = Vector{Int}(undef, nsites)
    cell_idx = Array{Int, D}(undef, sz...)

    for ci in CartesianIndices(sz)
        linear = LinearIndices(sz)[ci]
        cell_idx[ci] = (linear - 1) * K + 1
        for (α, s) in enumerate(sites)
            idx = (linear - 1) * K + α
            # position = sum of (ci[d]-1)*vec[d] + site position
            pos = ntuple(Val(D)) do d
                val = Float64(s[d])
                for dd in 1:D
                    val += (ci.I[dd] - 1) * Float64(vecs[dd][d])
                end
                val
            end
            site_positions[idx] = pos
            sublat[idx] = α
        end
    end

    # Build bonds — nearest neighbors based on lattice type
    bonds = _build_bonds(unitcell, sz, pbc, K, cell_idx)

    return FiniteLattice{D, typeof(unitcell)}(unitcell, sz, pbc, nsites, site_positions, bonds, cell_idx, sublat)
end

# Convenience constructors
FiniteLattice(unitcell::AbstractLattice{D}, sz::Vararg{Int, D}; kw...) where {D} =
    FiniteLattice(unitcell, sz; kw...)

"""
    site_index(fl::FiniteLattice, cell_coords::NTuple{D,Int}, sublattice::Int)

Get the global site index from unit cell grid coordinates (1-based) and sublattice index.
"""
function site_index(fl::FiniteLattice{D}, cell_coords::NTuple{D,Int}, sublattice::Int) where {D}
    K = nsites_per_cell(fl.unitcell)
    ci = CartesianIndex(cell_coords)
    linear = LinearIndices(fl.size)[ci]
    return (linear - 1) * K + sublattice
end

"""
    neighbors(fl::FiniteLattice, site::Int; label=nothing)

Return all neighbors of `site`. If `label` is specified, filter by bond label.
"""
function neighbors(fl::FiniteLattice, site::Int; label::Union{Nothing, Symbol} = nothing)
    result = Tuple{Int, Symbol}[]
    for b in fl.bonds
        if b.i == site
            (label === nothing || b.label == label) && push!(result, (b.j, b.label))
        elseif b.j == site
            (label === nothing || b.label == label) && push!(result, (b.i, b.label))
        end
    end
    return result
end

"""
    bonds_by_label(fl::FiniteLattice, label::Symbol)

Return all bonds of a given label as a vector of `(i, j)` pairs.
"""
function bonds_by_label(fl::FiniteLattice, label::Symbol)
    return [(b.i, b.j) for b in fl.bonds if b.label == label]
end

"""
    all_bond_labels(fl::FiniteLattice)

Return the set of all distinct bond labels in the lattice.
"""
all_bond_labels(fl::FiniteLattice) = unique(b.label for b in fl.bonds)

# ═══════════════════════════════════════════════════════════════
# Bond generation for specific lattice types
# ═══════════════════════════════════════════════════════════════

# Generic fallback: connect all nearest-neighbor sublattice pairs
# Subtype-specific methods below provide labeled bonds.

function _build_bonds(::AbstractLattice{D}, sz::NTuple{D,Int}, pbc::NTuple{D,Bool},
                      K::Int, cell_idx::Array{Int,D}) where {D}
    bonds = Bond[]
    # Default: just connect intra-cell sites with :nn label
    for ci in CartesianIndices(sz)
        base = cell_idx[ci]
        for α in 1:K, β in (α+1):K
            push!(bonds, Bond(base + α - 1, base + β - 1, :nn))
        end
    end
    return bonds
end

# ─────────── Honeycomb bonds (Kitaev-labeled: :x, :y, :z) ───────────
function _build_bonds(::HoneycombLattice, sz::NTuple{2,Int}, pbc::NTuple{2,Bool},
                      K::Int, cell_idx::Array{Int,2})
    Lx, Ly = sz
    bonds = Bond[]
    
    for ix in 1:Lx, iy in 1:Ly
        base = cell_idx[ix, iy]
        sA = base       # sublattice A (index 1)
        sB = base + 1   # sublattice B (index 2)
        
        # y-bond: intra-cell A-B
        push!(bonds, Bond(sA, sB, :y))

        # x-bond: B of (ix,iy) → A of (ix+1,iy)
        if ix < Lx
            base2 = cell_idx[ix + 1, iy]
            push!(bonds, Bond(sB, base2, :x))
        elseif pbc[1]
            base2 = cell_idx[1, iy]
            push!(bonds, Bond(sB, base2, :x))
        end

        # z-bond: B of (ix,iy) → A of (ix,iy+1)
        if iy < Ly
            base2 = cell_idx[ix, iy + 1]
            push!(bonds, Bond(sB, base2, :z))
        elseif pbc[2]
            base2 = cell_idx[ix, 1]
            push!(bonds, Bond(sA, base2 + 1, :z))  # A(ix,1) — B(ix,Ly) reversed: B(ix,iy)→A(ix,1)
            # Actually: z-bond with PBC wraps B(ix,Ly) → A(ix,1)
            # Correct: pop last and redo
            pop!(bonds)
            push!(bonds, Bond(sB, cell_idx[ix, 1], :z))
        end
    end
    return bonds
end

# ─────────── Square lattice bonds (:nn) ───────────
function _build_bonds(::SquareLattice, sz::NTuple{2,Int}, pbc::NTuple{2,Bool},
                      K::Int, cell_idx::Array{Int,2})
    Lx, Ly = sz
    bonds = Bond[]
    for ix in 1:Lx, iy in 1:Ly
        s = cell_idx[ix, iy]  # only 1 site per cell
        # +x direction
        if ix < Lx
            push!(bonds, Bond(s, cell_idx[ix + 1, iy], :nn))
        elseif pbc[1]
            push!(bonds, Bond(s, cell_idx[1, iy], :nn))
        end
        # +y direction  
        if iy < Ly
            push!(bonds, Bond(s, cell_idx[ix, iy + 1], :nn))
        elseif pbc[2]
            push!(bonds, Bond(s, cell_idx[ix, 1], :nn))
        end
    end
    return bonds
end

# ─────────── Triangular lattice bonds (:nn) ───────────
function _build_bonds(::TriangularLattice, sz::NTuple{2,Int}, pbc::NTuple{2,Bool},
                      K::Int, cell_idx::Array{Int,2})
    Lx, Ly = sz
    bonds = Bond[]
    _wrap(i, L, p) = p ? mod1(i, L) : (1 <= i <= L ? i : 0)
    
    for ix in 1:Lx, iy in 1:Ly
        s = cell_idx[ix, iy]
        # +a1 direction (ix+1, iy)
        jx = _wrap(ix + 1, Lx, pbc[1])
        jx > 0 && push!(bonds, Bond(s, cell_idx[jx, iy], :nn))
        # +a2 direction (ix, iy+1)
        jy = _wrap(iy + 1, Ly, pbc[2])
        jy > 0 && push!(bonds, Bond(s, cell_idx[ix, jy], :nn))
        # +a1-a2 direction (ix+1, iy-1) — the third NN for triangular
        jx2 = _wrap(ix + 1, Lx, pbc[1])
        jy2 = _wrap(iy - 1, Ly, pbc[2])
        (jx2 > 0 && jy2 > 0) && push!(bonds, Bond(s, cell_idx[jx2, jy2], :nn))
    end
    return bonds
end

# ─────────── Chain lattice bonds (:nn) ───────────
function _build_bonds(::ChainLattice, sz::Tuple{Int}, pbc::Tuple{Bool},
                      K::Int, cell_idx::Array{Int,1})
    L = sz[1]
    bonds = Bond[]
    for i in 1:L
        if i < L
            push!(bonds, Bond(cell_idx[i], cell_idx[i + 1], :nn))
        elseif pbc[1]
            push!(bonds, Bond(cell_idx[i], cell_idx[1], :nn))
        end
    end
    return bonds
end

# ─────────── Kagome lattice bonds (:nn) ───────────
function _build_bonds(::KagomeLattice, sz::NTuple{2,Int}, pbc::NTuple{2,Bool},
                      K::Int, cell_idx::Array{Int,2})
    Lx, Ly = sz
    bonds = Bond[]
    _wrap(i, L, p) = p ? mod1(i, L) : (1 <= i <= L ? i : 0)
    
    # K=3 sublattices: 1, 2, 3
    for ix in 1:Lx, iy in 1:Ly
        base = cell_idx[ix, iy]
        s1, s2, s3 = base, base + 1, base + 2
        
        # Intra-cell bonds
        push!(bonds, Bond(s1, s2, :nn))  # 1-2
        push!(bonds, Bond(s1, s3, :nn))  # 1-3
        
        # Inter-cell bonds
        # s2 connects to s3 of (ix+1, iy)
        jx = _wrap(ix + 1, Lx, pbc[1])
        jx > 0 && push!(bonds, Bond(s2, cell_idx[jx, iy] + 2, :nn))
        
        # s3 connects to s2 of (ix, iy+1)
        jy = _wrap(iy + 1, Ly, pbc[2])
        jy > 0 && push!(bonds, Bond(s3, cell_idx[ix, jy] + 1, :nn))
        
        # s3 connects to s2 of (ix+1, iy-1) — diagonal
        jx2 = _wrap(ix + 1, Lx, pbc[1])
        jy2 = _wrap(iy - 1, Ly, pbc[2])
        (jx2 > 0 && jy2 > 0) && push!(bonds, Bond(s3, cell_idx[jx2, jy2] + 1, :nn))
    end
    return bonds
end

# ═══════════════════════════════════════════════════════════════
# Translation symmetry permutations
# ═══════════════════════════════════════════════════════════════

"""
    translation_permutations(fl::FiniteLattice{D})

Return all non-identity translation permutations as `Vector{Vector{Int}}`.
Each inner vector `perm` satisfies: new_site[i] = old_site[perm[i]], i.e. site `i`
is mapped to site `perm[i]` under the translation.

For a lattice with `size = (L1, L2, ...)`, there are `prod(size) - 1` non-trivial translations.
"""
function translation_permutations(fl::FiniteLattice{D}) where {D}
    sz = fl.size
    K = nsites_per_cell(fl.unitcell)
    perms = Vector{Int}[]

    for shift in CartesianIndices(sz)
        # Skip identity (all shifts = 0)
        all(shift.I .== 1) && continue
        
        perm = Vector{Int}(undef, fl.nsites)
        for ci in CartesianIndices(sz)
            # Target cell after translation
            target = CartesianIndex(ntuple(d -> mod1(ci.I[d] + shift.I[d] - 1, sz[d]), Val(D)))
            src_base = fl.cell_indices[ci]
            dst_base = fl.cell_indices[target]
            for α in 1:K
                perm[src_base + α - 1] = dst_base + α - 1
            end
        end
        push!(perms, perm)
    end
    return perms
end

"""
    translation_group(fl::FiniteLattice{D})

Return the full translation group (including identity) as `Vector{Vector{Int}}`.
"""
function translation_group(fl::FiniteLattice{D}) where {D}
    K = nsites_per_cell(fl.unitcell)
    identity = collect(1:fl.nsites)
    return vcat([identity], translation_permutations(fl))
end
