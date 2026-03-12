# ═══════════════════════════════════════════════════════════════
# PointGroupSymmetry: rotations, reflections
# (Placeholder for Phase 4 — not yet implemented)
# ═══════════════════════════════════════════════════════════════

"""
    PointGroupSymmetry <: AbstractSpatialSymmetry

Point group symmetry (rotations, reflections) of the finite lattice.
To be implemented in Phase 4.

# Fields
- `permutations::Vector{Vector{Int}}` — point group operation permutations
- `characters::Vector{ComplexF64}` — characters of the target irrep
"""
struct PointGroupSymmetry <: AbstractSpatialSymmetry
    permutations::Vector{Vector{Int}}
    characters::Vector{ComplexF64}
end
