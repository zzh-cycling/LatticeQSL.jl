# ═══════════════════════════════════════════════════════════════
# Predefined lattice types
# ═══════════════════════════════════════════════════════════════

# ─────────── Honeycomb ───────────
"""
    struct HoneycombLattice <: AbstractLattice{2}

Type representing 2D Honeycomb Lattice (2 sites per unit cell).
"""
struct HoneycombLattice <: AbstractLattice{2} end

lattice_vectors(::HoneycombLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))
reciprocal_lattice_vector(::HoneycombLattice) = ((2π, -2π / sqrt(3)), (0.0, 4π / sqrt(3)))
lattice_sites(::HoneycombLattice) = ((0.0, 0.0), (0.5, 0.5 / sqrt(3)))
reciprocal_lattice_sites(::HoneycombLattice) = ((0.0, 0.0), (1 / 3, 1 / sqrt(3)))

# ─────────── Square ───────────
"""
    struct SquareLattice <: AbstractLattice{2}

Type representing 2D Square Lattice (1 site per unit cell).
"""
struct SquareLattice <: AbstractLattice{2} end

lattice_vectors(::SquareLattice) = ((1.0, 0.0), (0.0, 1.0))
reciprocal_lattice_vector(::SquareLattice) = ((2π, 0.0), (0.0, 2π))
lattice_sites(::SquareLattice) = ((0.0, 0.0),)
reciprocal_lattice_sites(::SquareLattice) = ((0.0, 0.0),)

# ─────────── Triangular ───────────
"""
    struct TriangularLattice <: AbstractLattice{2}

Type representing 2D Triangular Lattice (1 site per unit cell).
"""
struct TriangularLattice <: AbstractLattice{2} end

lattice_vectors(::TriangularLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))
reciprocal_lattice_vector(::TriangularLattice) = ((2π, -2π / sqrt(3)), (0.0, 4π / sqrt(3)))
lattice_sites(::TriangularLattice) = ((0.0, 0.0),)
reciprocal_lattice_sites(::TriangularLattice) = ((0.0, 0.0),)

# ─────────── Chain ───────────
"""
    struct ChainLattice <: AbstractLattice{1}

Type representing 1D Chain Lattice (1 site per unit cell).
"""
struct ChainLattice <: AbstractLattice{1} end

lattice_vectors(::ChainLattice) = ((1.0,),)
reciprocal_lattice_vector(::ChainLattice) = ((2π,),)
lattice_sites(::ChainLattice) = ((0.0,),)
reciprocal_lattice_sites(::ChainLattice) = ((0.0,),)

# ─────────── Kagome ───────────
"""
    struct KagomeLattice <: AbstractLattice{2}

Type representing 2D Kagome Lattice (3 sites per unit cell).
"""
struct KagomeLattice <: AbstractLattice{2} end

lattice_vectors(::KagomeLattice) = ((1.0, 0.0), (0.5, 0.5 * sqrt(3)))
reciprocal_lattice_vector(::KagomeLattice) = ((2π, -2π / sqrt(3)), (0.0, 4π / sqrt(3)))
lattice_sites(::KagomeLattice) = ((0.0, 0.0), (0.25, 0.25 * sqrt(3)), (0.75, 0.25 * sqrt(3)))
reciprocal_lattice_sites(::KagomeLattice) = ((0.0, 0.0), (1 / 3, 1 / sqrt(3)), (2 / 3, 1 / sqrt(3)))

# ─────────── Lieb ───────────
"""
    struct LiebLattice <: AbstractLattice{2}

Type representing 2D Lieb Lattice (3 sites per unit cell).
"""
struct LiebLattice <: AbstractLattice{2} end

lattice_vectors(::LiebLattice) = ((1.0, 0.0), (0.0, 1.0))
reciprocal_lattice_vector(::LiebLattice) = ((2π, 0.0), (0.0, 2π))
lattice_sites(::LiebLattice) = ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5))

# ─────────── Rectangular ───────────
"""
    struct RectangularLattice <: AbstractLattice{2}

Type representing 2D Rectangular Lattice with adjustable aspect ratio (1 site per unit cell).
"""
struct RectangularLattice <: AbstractLattice{2}
    aspect_ratio::Float64
end

lattice_vectors(r::RectangularLattice) = ((1.0, 0.0), (0.0, r.aspect_ratio))
reciprocal_lattice_vector(r::RectangularLattice) = ((2π, 0.0), (0.0, 2π / r.aspect_ratio))
lattice_sites(::RectangularLattice) = ((0.0, 0.0),)
reciprocal_lattice_sites(::RectangularLattice) = ((0.0, 0.0),)

# ─────────── GeneralLattice ───────────
"""
    GeneralLattice{D}(vectors, sites)

A user-defined lattice with arbitrary vectors and sites.
"""
struct GeneralLattice{D} <: AbstractLattice{D}
    vectors::NTuple{D, NTuple{D, Float64}}
    sites::Vector{NTuple{D, Float64}}
    function GeneralLattice(vectors, sites)
        D = length(vectors)
        vecs = ntuple(i -> Tuple(Float64.(vectors[i])), Val(D))
        ss = [Tuple(Float64.(s)) for s in sites]
        new{D}(vecs, ss)
    end
end

lattice_vectors(g::GeneralLattice) = g.vectors
lattice_sites(g::GeneralLattice) = Tuple(g.sites)
function reciprocal_lattice_vector(g::GeneralLattice{D}) where {D}
    M = reduce(hcat, [collect(v) for v in g.vectors])
    invM = inv(M)'
    return ntuple(i -> Tuple(2π .* invM[:, i]), Val(D))
end
function reciprocal_lattice_sites(g::GeneralLattice{D}) where {D}
    M = reduce(hcat, [collect(v) for v in g.vectors])
    invM = inv(M)
    return ntuple(k -> Tuple(2π .* (invM * collect(g.sites[k]))), Val(length(g.sites)))
end
