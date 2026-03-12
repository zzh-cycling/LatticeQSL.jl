# ═══════════════════════════════════════════════════════════════
# AbstractModel: interface for spin models
# ═══════════════════════════════════════════════════════════════

"""
    AbstractModel

Abstract supertype for spin Hamiltonian models.

All models must implement:
- `model_bonds(model, fl::FiniteLattice)` → Vector of (i, j, coeff, bond_type) 
- `apply_model_bond(model, state, i, j, bond_type)` → Vector{(new_state, coeff)}
"""
abstract type AbstractModel end

"""
    HeisenbergModel <: AbstractModel

Isotropic Heisenberg model: H = J Σ_{⟨i,j⟩} Sᵢ·Sⱼ

# Fields
- `J::Float64` — coupling constant (J > 0 = antiferromagnetic, J < 0 = ferromagnetic)
- `bond_label::Symbol` — which bonds to include (default `:nn`)
"""
struct HeisenbergModel <: AbstractModel
    J::Float64
    bond_label::Symbol
end
HeisenbergModel(; J::Real = 1.0, bond_label::Symbol = :nn) = HeisenbergModel(Float64(J), bond_label)

"""
    KitaevModel <: AbstractModel

Kitaev model: H = Σ_{⟨i,j⟩_α} Jα Sᵢᵅ Sⱼᵅ, where α ∈ {x, y, z} depends on bond type.

# Fields
- `Jx::Float64`, `Jy::Float64`, `Jz::Float64` — coupling constants
"""
struct KitaevModel <: AbstractModel
    Jx::Float64
    Jy::Float64
    Jz::Float64
end
KitaevModel(; Jx::Real = 1.0, Jy::Real = 1.0, Jz::Real = 1.0) =
    KitaevModel(Float64(Jx), Float64(Jy), Float64(Jz))

"""
    IsingModel <: AbstractModel

Transverse-field Ising model: H = J Σ_{⟨i,j⟩} SᶻᵢSᶻⱼ + h Σᵢ Sˣᵢ

# Fields
- `J::Float64` — Ising coupling
- `h::Float64` — transverse field
- `bond_label::Symbol` — which bonds to include (default `:nn`)
"""
struct IsingModel <: AbstractModel
    J::Float64
    h::Float64
    bond_label::Symbol
end
IsingModel(; J::Real = 1.0, h::Real = 0.0, bond_label::Symbol = :nn) =
    IsingModel(Float64(J), Float64(h), bond_label)

"""
    XXZModel <: AbstractModel

XXZ model: H = Σ_{⟨i,j⟩} [Jxy(SˣᵢSˣⱼ + SʸᵢSʸⱼ) + Jz SᶻᵢSᶻⱼ]

# Fields
- `Jxy::Float64` — XX and YY coupling
- `Jz::Float64` — ZZ coupling  
- `bond_label::Symbol`
"""
struct XXZModel <: AbstractModel
    Jxy::Float64
    Jz::Float64
    bond_label::Symbol
end
XXZModel(; Jxy::Real = 1.0, Jz::Real = 1.0, bond_label::Symbol = :nn) =
    XXZModel(Float64(Jxy), Float64(Jz), bond_label)

"""
    GeneralSpinModel <: AbstractModel

General spin model: H = Σ_{⟨i,j⟩} [Jx SˣᵢSˣⱼ + Jy SʸᵢSʸⱼ + Jz SᶻᵢSᶻⱼ]

Can represent Heisenberg (Jx=Jy=Jz), XXZ (Jx=Jy≠Jz), XY (Jz=0), Ising (Jx=Jy=0), etc.

# Fields
- `Jx::Float64`, `Jy::Float64`, `Jz::Float64` — coupling constants for each component
- `bond_label::Symbol` — which bonds to include
"""
struct GeneralSpinModel <: AbstractModel
    Jx::Float64
    Jy::Float64
    Jz::Float64
    bond_label::Symbol
end
GeneralSpinModel(; Jx::Real = 1.0, Jy::Real = 1.0, Jz::Real = 1.0, bond_label::Symbol = :nn) =
    GeneralSpinModel(Float64(Jx), Float64(Jy), Float64(Jz), bond_label)
