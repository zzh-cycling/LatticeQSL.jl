# ═══════════════════════════════════════════════════════════════
# Model implementations: Heisenberg, Kitaev, Ising, GeneralSpin
# ═══════════════════════════════════════════════════════════════

"""
    model_bonds(model::AbstractModel, fl::FiniteLattice)

Return all bonds contributing to the Hamiltonian as a vector of
`(i, j, coupling, bond_type)` tuples.
"""
function model_bonds end

"""
    apply_model_bond(model::AbstractModel, state, i, j, bond_type)

Apply one bond term to a basis state.
Returns a vector of `(new_state, coefficient)` pairs.
"""
function apply_model_bond end

# ─────────── Heisenberg ───────────
function model_bonds(model::HeisenbergModel, fl::FiniteLattice)
    return [(i, j, model.J, model.bond_label) for (i, j) in bonds_by_label(fl, model.bond_label)]
end

function apply_model_bond(::HeisenbergModel, state::BitStr{N, Int}, i::Int, j::Int, bond_type::Symbol) where {N}
    return apply_heisenberg(state, i, j)
end

# ─────────── Kitaev ───────────
function model_bonds(model::KitaevModel, fl::FiniteLattice)
    bonds = Tuple{Int, Int, Float64, Symbol}[]
    for b in fl.bonds
        coupling = if b.label == :x
            model.Jx
        elseif b.label == :y
            model.Jy
        elseif b.label == :z
            model.Jz
        else
            continue  # skip non-Kitaev bonds
        end
        push!(bonds, (b.i, b.j, coupling, b.label))
    end
    return bonds
end

function apply_model_bond(::KitaevModel, state::BitStr{N, Int}, i::Int, j::Int, bond_type::Symbol) where {N}
    return apply_kitaev_bond(state, i, j, bond_type)
end

# ─────────── Ising ───────────
function model_bonds(model::IsingModel, fl::FiniteLattice)
    # IsingModel uses bonds + a transverse field handled separately in the builder
    return [(i, j, model.J, model.bond_label) for (i, j) in bonds_by_label(fl, model.bond_label)]
end

function apply_model_bond(::IsingModel, state::BitStr{N, Int}, i::Int, j::Int, bond_type::Symbol) where {N}
    _, coeff = apply_SzSz(state, i, j)
    return [(state, coeff)]
end

# ─────────── XXZ ───────────
function model_bonds(model::XXZModel, fl::FiniteLattice)
    return [(i, j, 1.0, model.bond_label) for (i, j) in bonds_by_label(fl, model.bond_label)]
end

function apply_model_bond(model::XXZModel, state::BitStr{N, Int}, i::Int, j::Int, bond_type::Symbol) where {N}
    results = Tuple{BitStr{N, Int}, ComplexF64}[]
    
    # Jz SzSz
    _, coeff_zz = apply_SzSz(state, i, j)
    coeff_zz != 0 && push!(results, (state, model.Jz * coeff_zz))
    
    # Jxy/2 (S⁺ᵢS⁻ⱼ + S⁻ᵢS⁺ⱼ)
    r1 = apply_SpSm(state, i, j)
    r1 !== nothing && push!(results, (r1[1], model.Jxy * 0.5 * r1[2]))
    
    r2 = apply_SmSp(state, i, j)
    r2 !== nothing && push!(results, (r2[1], model.Jxy * 0.5 * r2[2]))
    
    return results
end

# ─────────── GeneralSpin ───────────
function model_bonds(model::GeneralSpinModel, fl::FiniteLattice)
    return [(i, j, 1.0, model.bond_label) for (i, j) in bonds_by_label(fl, model.bond_label)]
end

function apply_model_bond(model::GeneralSpinModel, state::BitStr{N, Int}, i::Int, j::Int, bond_type::Symbol) where {N}
    results = Tuple{BitStr{N, Int}, ComplexF64}[]
    
    # Jz SzSz
    _, coeff_zz = apply_SzSz(state, i, j)
    coeff_zz != 0 && push!(results, (state, model.Jz * coeff_zz))
    
    # Jx SxSx = Jx * (1/4)σxσx: flip both spins
    if model.Jx != 0
        new_state = flip_spin(flip_spin(state, i), j)
        push!(results, (new_state, ComplexF64(model.Jx * 0.25)))
    end
    
    # Jy SySy = Jy * (1/4)σyσy
    if model.Jy != 0
        si = Int((state.buf >> (i - 1)) & 0x1)
        sj = Int((state.buf >> (j - 1)) & 0x1)
        new_state = flip_spin(flip_spin(state, i), j)
        ci = si == 0 ? im : -im
        cj = sj == 0 ? im : -im
        coeff = model.Jy * 0.25 * ci * cj
        push!(results, (new_state, ComplexF64(coeff)))
    end
    
    return results
end
