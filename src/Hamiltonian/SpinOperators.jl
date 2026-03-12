# ═══════════════════════════════════════════════════════════════
# SpinOperators: efficient spin-1/2 operator actions on BitStr states
#
# All operators act on a single state and return (new_state, coefficient) pairs.
# Site indices are 1-based.
# Convention: Sᵢ = (1/2)σᵢ where σ are Pauli matrices.
# ═══════════════════════════════════════════════════════════════

# ─────────── Single-site operators ───────────

"""
    apply_Sz(state::BitStr{N}, site::Int) -> (state, coeff)

Apply Sᶻᵢ = (1/2)σᶻᵢ. Diagonal: eigenvalue = +1/2 (up) or -1/2 (down).
"""
function apply_Sz(state::BitStr{N, Int}, site::Int) where {N}
    bit = Int((state.buf >> (site - 1)) & 0x1)
    return (state, 0.5 * (2 * bit - 1))
end

"""
    apply_Sp(state::BitStr{N}, site::Int) -> (new_state, coeff) or nothing

Apply S⁺ᵢ = Sˣᵢ + iSʸᵢ. Returns nothing if spin is already up.
S⁺|↓⟩ = |↑⟩, S⁺|↑⟩ = 0.
"""
function apply_Sp(state::BitStr{N, Int}, site::Int) where {N}
    if ((state.buf >> (site - 1)) & 0x1) == 0  # spin down → flip to up
        return (flip_spin(state, site), 1.0)
    end
    return nothing  # already up, S⁺|↑⟩ = 0
end

"""
    apply_Sm(state::BitStr{N}, site::Int) -> (new_state, coeff) or nothing

Apply S⁻ᵢ = Sˣᵢ - iSʸᵢ. Returns nothing if spin is already down.
S⁻|↑⟩ = |↓⟩, S⁻|↓⟩ = 0.
"""
function apply_Sm(state::BitStr{N, Int}, site::Int) where {N}
    if ((state.buf >> (site - 1)) & 0x1) == 1  # spin up → flip to down
        return (flip_spin(state, site), 1.0)
    end
    return nothing  # already down, S⁻|↓⟩ = 0
end

"""
    apply_Sx(state::BitStr{N}, site::Int) -> (new_state, coeff)

Apply Sˣᵢ = (1/2)(S⁺ + S⁻) = (1/2)σˣᵢ. Always flips the spin.
"""
function apply_Sx(state::BitStr{N, Int}, site::Int) where {N}
    return (flip_spin(state, site), 0.5)
end

"""
    apply_Sy(state::BitStr{N}, site::Int) -> (new_state, coeff)

Apply Sʸᵢ = (1/2i)(S⁺ - S⁻) = (1/2)σʸᵢ.
σʸ|↓⟩ = i|↑⟩, σʸ|↑⟩ = -i|↓⟩
So Sʸ|↓⟩ = (i/2)|↑⟩, Sʸ|↑⟩ = (-i/2)|↓⟩
"""
function apply_Sy(state::BitStr{N, Int}, site::Int) where {N}
    if readbit(state, site) == 0  # down
        return (flip_spin(state, site), im * 0.5)
    else  # up
        return (flip_spin(state, site), -im * 0.5)
    end
end

# ─────────── Two-site operators ───────────

"""
    apply_SzSz(state::BitStr{N}, i::Int, j::Int) -> (state, coeff)

Apply SᶻᵢSᶻⱼ. Diagonal operator.
"""
function apply_SzSz(state::BitStr{N, Int}, i::Int, j::Int) where {N}
    si = 2 * Int((state.buf >> (i - 1)) & 0x1) - 1  # ±1
    sj = 2 * Int((state.buf >> (j - 1)) & 0x1) - 1  # ±1
    return (state, 0.25 * si * sj)
end

"""
    apply_SpSm(state::BitStr{N}, i::Int, j::Int) -> (new_state, coeff) or nothing

Apply S⁺ᵢS⁻ⱼ. Flips spin i up and spin j down. Returns nothing if not possible.
"""
function apply_SpSm(state::BitStr{N, Int}, i::Int, j::Int) where {N}
    # Need: site i is down (to flip up), site j is up (to flip down)
    if ((state.buf >> (i - 1)) & 0x1) == 0 && ((state.buf >> (j - 1)) & 0x1) == 1
        new_state = flip_spin(flip_spin(state, i), j)
        return (new_state, 1.0)
    end
    return nothing
end

"""
    apply_SmSp(state::BitStr{N}, i::Int, j::Int) -> (new_state, coeff) or nothing

Apply S⁻ᵢS⁺ⱼ. Flips spin i down and spin j up. Returns nothing if not possible.
"""
function apply_SmSp(state::BitStr{N, Int}, i::Int, j::Int) where {N}
    if ((state.buf >> (i - 1)) & 0x1) == 1 && ((state.buf >> (j - 1)) & 0x1) == 0
        new_state = flip_spin(flip_spin(state, i), j)
        return (new_state, 1.0)
    end
    return nothing
end

"""
    apply_heisenberg(state::BitStr{N}, i::Int, j::Int) -> Vector{Tuple{BitStr{N}, Float64}}

Apply Sᵢ·Sⱼ = SᶻᵢSᶻⱼ + (1/2)(S⁺ᵢS⁻ⱼ + S⁻ᵢS⁺ⱼ) to a state.
Returns a list of (new_state, coefficient) pairs.
"""
function apply_heisenberg(state::BitStr{N, Int}, i::Int, j::Int) where {N}
    results = Tuple{BitStr{N, Int}, Float64}[]
    
    # SzSz (diagonal)
    _, coeff_zz = apply_SzSz(state, i, j)
    coeff_zz != 0 && push!(results, (state, coeff_zz))
    
    # (1/2) S⁺ᵢS⁻ⱼ
    r1 = apply_SpSm(state, i, j)
    r1 !== nothing && push!(results, (r1[1], 0.5 * r1[2]))
    
    # (1/2) S⁻ᵢS⁺ⱼ
    r2 = apply_SmSp(state, i, j)
    r2 !== nothing && push!(results, (r2[1], 0.5 * r2[2]))
    
    return results
end

"""
    apply_kitaev_bond(state::BitStr{N}, i::Int, j::Int, bond_type::Symbol) -> Vector{Tuple{BitStr{N}, ComplexF64}}

Apply the Kitaev interaction SᵅᵢSᵅⱼ for α ∈ {:x, :y, :z} to a state.
Returns a list of (new_state, coefficient) pairs.
"""
function apply_kitaev_bond(state::BitStr{N, Int}, i::Int, j::Int, bond_type::Symbol) where {N}
    if bond_type == :z
        _, coeff = apply_SzSz(state, i, j)
        return [(state, ComplexF64(coeff))]
    elseif bond_type == :x
        # SxSx = (1/4)σxσx: flip both spins, coefficient = 1/4
        new_state = flip_spin(flip_spin(state, i), j)
        return [(new_state, ComplexF64(0.25))]
    elseif bond_type == :y
        # SySy = (1/4)σyσy
        # σy|0⟩ = i|1⟩, σy|1⟩ = -i|0⟩
        # σy_i σy_j |si sj⟩ → coefficient depends on si, sj
        si = Int((state.buf >> (i - 1)) & 0x1)
        sj = Int((state.buf >> (j - 1)) & 0x1)
        new_state = flip_spin(flip_spin(state, i), j)
        # Phase: (i)^(1-si) × (-i)^si × (i)^(1-sj) × (-i)^sj
        # = i^(1-2si) × i^(1-2sj) = i^(2 - 2si - 2sj)  ... actually:
        # σy|0⟩ = i|1⟩  → coeff = i
        # σy|1⟩ = -i|0⟩ → coeff = -i
        ci = si == 0 ? im : -im
        cj = sj == 0 ? im : -im
        coeff = 0.25 * ci * cj  # (1/2)σy × (1/2)σy = (1/4)σyσy
        return [(new_state, ComplexF64(coeff))]
    else
        error("Unknown Kitaev bond type: $bond_type. Expected :x, :y, or :z")
    end
end
