# ═══════════════════════════════════════════════════════════════
# Symmetry: abstract interface for symmetry operations
# ═══════════════════════════════════════════════════════════════

"""
    AbstractSymmetry

Abstract type for symmetry operations that can be used to reduce the Hilbert space.
"""
abstract type AbstractSymmetry end

"""
    AbstractSpatialSymmetry <: AbstractSymmetry

Spatial symmetries (translations, rotations, reflections) that act as permutations on sites.
"""
abstract type AbstractSpatialSymmetry <: AbstractSymmetry end

"""
    AbstractInternalSymmetry <: AbstractSymmetry

Internal symmetries (Sz conservation, spin inversion) that act on the spin degrees of freedom.
"""
abstract type AbstractInternalSymmetry <: AbstractSymmetry end
