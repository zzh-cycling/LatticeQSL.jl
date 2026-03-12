# ═══════════════════════════════════════════════════════════════
# TranslationSymmetry: momentum-space basis construction
# ═══════════════════════════════════════════════════════════════

"""
    TranslationSymmetry <: AbstractSpatialSymmetry

Translation symmetry group of a finite lattice. By constructing a symmetry-adapted basis,
the Hilbert space is decomposed into momentum sectors labeled by **k**.

# Fields
- `permutations::Vector{Vector{Int}}` — all translation permutations (including identity)
- `momentum_index::Int` — index of the target momentum sector (1-based, in the order of `CartesianIndices(lattice.size)`)
"""
struct TranslationSymmetry <: AbstractSpatialSymmetry
    permutations::Vector{Vector{Int}}
    momentum_index::Int
end

"""
    TranslationSymmetry(fl::FiniteLattice, k_index::Int)

Create a TranslationSymmetry from a FiniteLattice and a momentum sector index.
"""
function TranslationSymmetry(fl::FiniteLattice, k_index::Int)
    perms = translation_group(fl)
    return TranslationSymmetry(perms, k_index)
end

"""
    momentum_values(fl::FiniteLattice{D}) -> Vector{NTuple{D, Float64}}

Return all allowed momentum values k = (2π n₁/L₁, 2π n₂/L₂, ...) for the finite lattice.
The ordering matches `CartesianIndices(fl.size)`.
"""
function momentum_values(fl::FiniteLattice{D}) where {D}
    sz = fl.size
    kvalues = NTuple{D, Float64}[]
    for ci in CartesianIndices(sz)
        k = ntuple(d -> 2π * (ci.I[d] - 1) / sz[d], Val(D))
        push!(kvalues, k)
    end
    return kvalues
end
