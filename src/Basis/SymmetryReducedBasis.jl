# ═══════════════════════════════════════════════════════════════
# SymmetryReducedBasis: construct symmetry-adapted basis states
#
# Given a symmetry group G with elements {g₁, g₂, ...} and an irrep χ,
# the symmetry-adapted state for representative |α⟩ is:
#
#   |k, α⟩ = (1/𝒩_α) Σ_g χ*(g) g|α⟩
#
# where 𝒩_α is a normalization factor.
#
# We store:
# - representative states
# - normalization factors
# - a mapping from any state to its representative's index
# ═══════════════════════════════════════════════════════════════

"""
    SymmetryReducedBasis{N}

A symmetry-reduced basis for N spin-1/2 sites.

# Fields
- `nsites::Int` — number of sites
- `parent_basis::SpinBasis{N}` — the original (possibly Sz-conserved) basis
- `representatives::Vector{BitStr{N, Int}}` — representative states (sorted)
- `rep_indices::Vector{Int}` — index into `representatives` for each state in parent_basis (0 if not a representative)
- `norms::Vector{Float64}` — normalization factors for each representative
- `orbit_sizes::Vector{Int}` — number of unique states in each orbit
- `dim::Int` — dimension of the reduced basis
- `characters::Vector{ComplexF64}` — characters χ(g) for each group element
- `permutations::Vector{Vector{Int}}` — group elements as site permutations
"""
struct SymmetryReducedBasis{N}
    nsites::Int
    parent_basis::SpinBasis{N}
    representatives::Vector{BitStr{N, Int}}
    state_to_rep::Dict{BitStr{N, Int}, Int}  # maps any state → index of its representative in `representatives`
    state_phases::Dict{BitStr{N, Int}, ComplexF64}  # phase relating a state to its orbit representative
    norms::Vector{Float64}  # √(effective orbit weight)
    orbit_sizes::Vector{Int}
    dim::Int
    characters::Vector{ComplexF64}
    permutations::Vector{Vector{Int}}
end

"""
    SymmetryReducedBasis(parent_basis::SpinBasis{N}, sym::TranslationSymmetry, fl::FiniteLattice)

Construct a translation-symmetry-adapted basis for a given momentum sector.

The characters are χ(Tₙ) = exp(-i k·Rₙ) where Rₙ is the translation vector.
"""
function SymmetryReducedBasis(parent_basis::SpinBasis{N}, sym::TranslationSymmetry,
                              fl::FiniteLattice{D}) where {N, D}
    perms = sym.permutations
    ngroup = length(perms)
    
    # Compute characters: χ(g) = exp(-i k · R_g)
    kvalues = momentum_values(fl)
    k = kvalues[sym.momentum_index]
    
    # Characters for each group element
    # The group elements are indexed by CartesianIndices(fl.size)
    characters = Vector{ComplexF64}(undef, ngroup)
    sz = fl.size
    idx = 0
    for ci in CartesianIndices(sz)
        idx += 1
        phase = sum(k[d] * (ci.I[d] - 1) for d in 1:D)
        characters[idx] = exp(im * phase)
    end

    return _build_translation_reduced_basis(parent_basis, perms, characters)
end

"""
    SymmetryReducedBasis(parent_basis::SpinBasis{N}, sym::TranslationSymmetry,
                         fl::FiniteLattice{1}, inversion::Int)

Construct a 1D translation + reflection reduced basis. This follows the same
sector-selection logic as the legacy `iso_K2MSS` implementation: only the
`k = 0` and `k = π` sectors are reflection-resolved.
"""
function SymmetryReducedBasis(parent_basis::SpinBasis{N}, sym::TranslationSymmetry,
                              fl::FiniteLattice{1}, inversion::Int) where {N}
    @assert inversion == 1 || inversion == -1 "inversion must be ±1, got $inversion"

    k = sym.momentum_index - 1
    @assert k == 0 || k == fl.size[1] ÷ 2 "reflection reduction requires k = 0 or π sector"

    tbasis = SymmetryReducedBasis(parent_basis, sym, fl)
    return _build_reflection_reduced_basis(tbasis, inversion, fl.size[1])
end

"""
    SymmetryReducedBasis(parent_basis::SpinBasis{N}, perms::Vector{Vector{Int}}, characters::Vector{ComplexF64})

General constructor: build a symmetry-reduced basis from explicit permutations and characters.
"""
function SymmetryReducedBasis(parent_basis::SpinBasis{N}, perms::Vector{Vector{Int}},
                              characters::Vector{ComplexF64}) where {N}
    return _build_reduced_basis(parent_basis, perms, characters)
end


function _build_reduced_basis(parent_basis::SpinBasis{N}, perms::Vector{Vector{Int}},
                              characters::Vector{ComplexF64}) where {N}
    T = BitStr{N, Int}

    visited = Set{T}()
    representatives = T[]
    state_to_rep = Dict{T, Int}()
    state_phases = Dict{T, ComplexF64}()
    norms = Float64[]
    orbit_sizes = Int[]

    for state in parent_basis.states

        state in visited && continue

        orbit_map = Dict{T, Int}()
        for (ig, perm) in enumerate(perms)
            orbit_map[permute_bits(state, perm)] = get(orbit_map, permute_bits(state, perm), ig)
        end
        orbit_states = sort!(collect(keys(orbit_map)))
        union!(visited, orbit_states)

        rep = first(orbit_states)
        push!(representatives, rep)
        push!(orbit_sizes, length(orbit_states))
        push!(norms, sqrt(length(orbit_states)))
        rep_idx = length(representatives)

        for orbit_state in orbit_states
            state_to_rep[orbit_state] = rep_idx
            state_phases[orbit_state] = characters[_perm_index_to_representative(orbit_state, rep, perms)]
        end
    end

    dim = length(representatives)

    return SymmetryReducedBasis{N}(
        parent_basis.nsites,
        parent_basis,
        representatives,
        state_to_rep,
        state_phases,
        norms,
        orbit_sizes,
        dim,
        characters,
        perms
    )
end

function _build_translation_reduced_basis(parent_basis::SpinBasis{N}, perms::Vector{Vector{Int}},
                                          characters::Vector{ComplexF64}) where {N}
    T = BitStr{N, Int}
    visited = Set{T}()
    representatives = T[]
    state_to_rep = Dict{T, Int}()
    state_phases = Dict{T, ComplexF64}()
    norms = Float64[]
    orbit_sizes = Int[]

    for state in parent_basis.states
        state in visited && continue

        orbit_states = _orbit_states(state, perms)
        union!(visited, orbit_states)
        rep = first(orbit_states)

        stabilizer_phase = zero(ComplexF64)
        for (ig, perm) in enumerate(perms)
            if permute_bits(rep, perm) == rep
                stabilizer_phase += characters[ig]
            end
        end

        abs(stabilizer_phase) < 1e-12 && continue

        push!(representatives, rep)
        push!(orbit_sizes, length(orbit_states))
        push!(norms, sqrt(length(orbit_states)))
        rep_idx = length(representatives)

        for orbit_state in orbit_states
            state_to_rep[orbit_state] = rep_idx
            state_phases[orbit_state] = characters[_perm_index_to_representative(orbit_state, rep, perms)]
        end
    end

    return SymmetryReducedBasis{N}(
        parent_basis.nsites,
        parent_basis,
        representatives,
        state_to_rep,
        state_phases,
        norms,
        orbit_sizes,
        length(representatives),
        characters,
        perms,
    )
end

function _build_reflection_reduced_basis(tbasis::SymmetryReducedBasis{N}, inversion::Int,
                                         ntranslations::Int) where {N}
    T = BitStr{N, Int}
    k = _translation_sector_from_characters(tbasis.characters, ntranslations)

    representatives = T[]
    state_to_rep = Dict{T, Int}()
    state_phases = Dict{T, ComplexF64}()
    norms = Float64[]
    orbit_sizes = Int[]
    seen_reps = Set{T}()

    keep_fixed = (inversion == 1 && k == 0) || (inversion == -1 && k == div(ntranslations, 2))

    for (idx, rep) in enumerate(tbasis.representatives)
        rep in seen_reps && continue

        reflected = _translation_representative(BitBasis.breflect(rep), tbasis.permutations)
        canonical = min(rep, reflected)

        if keep_fixed
            rep == canonical || continue
            q = rep == reflected ? 1 : 2
        else
            rep == reflected && continue
            rep == canonical || continue
            q = 2
        end

        push!(representatives, canonical)
        push!(orbit_sizes, q * tbasis.orbit_sizes[idx])
        push!(norms, sqrt(q * tbasis.orbit_sizes[idx]))
        rep_idx = length(representatives)

        for (state, state_rep_idx) in tbasis.state_to_rep
            krep = tbasis.representatives[state_rep_idx]
            if krep == rep || krep == reflected
                state_to_rep[state] = rep_idx
                state_phases[state] = tbasis.state_phases[state]
            end
        end

        push!(seen_reps, rep)
        push!(seen_reps, reflected)
    end

    return SymmetryReducedBasis{N}(
        tbasis.nsites,
        tbasis.parent_basis,
        representatives,
        state_to_rep,
        state_phases,
        norms,
        orbit_sizes,
        length(representatives),
        tbasis.characters,
        tbasis.permutations,
    )
end

function _orbit_states(state::BitStr{N, Int}, perms::Vector{Vector{Int}}) where {N}
    orbit = Set{BitStr{N, Int}}()
    for perm in perms
        push!(orbit, permute_bits(state, perm))
    end
    return sort!(collect(orbit))
end

function _translation_representative(state::BitStr{N, Int}, perms::Vector{Vector{Int}}) where {N}
    return first(_orbit_states(state, perms))
end

function _perm_index_to_representative(state::BitStr{N, Int}, rep::BitStr{N, Int},
                                       perms::Vector{Vector{Int}}) where {N}
    for (ig, perm) in enumerate(perms)
        if permute_bits(state, perm) == rep
            return ig
        end
    end
    error("Failed to map state to its representative")
end

function _translation_sector_from_characters(characters::Vector{ComplexF64}, ntranslations::Int)
    phase = angle(characters[min(2, length(characters))])
    k = round(Int, phase * ntranslations / (2π))
    return mod(k, ntranslations)
end

"""
    representative_index(basis::SymmetryReducedBasis{N}, state::BitStr{N, Int}) -> Int

Find the index of the representative of `state`'s orbit. Returns 0 if not in this sector.
"""
function representative_index(basis::SymmetryReducedBasis{N}, state::BitStr{N, Int}) where {N}
    return get(basis.state_to_rep, state, 0)
end

"""
    matrix_element_factor(basis::SymmetryReducedBasis{N}, 
                          rep_a::BitStr{N,Int}, idx_a::Int,
                          new_state::BitStr{N,Int}) -> (Int, ComplexF64)

Given that the Hamiltonian maps representative |a⟩ to some state |b'⟩ = new_state,
find the representative index of |b'⟩ and the phase factor needed for the matrix element.

Returns `(idx_b, factor)` where:
- `idx_b` is the index of the representative of new_state (0 if outside sector)
- `factor` includes the character phase and normalization: χ(g) × 𝒩_a / 𝒩_b
  where g maps the representative of b's orbit to new_state.
"""
function matrix_element_factor(basis::SymmetryReducedBasis{N},
                               idx_a::Int, new_state::BitStr{N, Int}) where {N}
    idx_b = representative_index(basis, new_state)
    idx_b == 0 && return (0, zero(ComplexF64))

    phase = get(basis.state_phases, new_state, zero(ComplexF64))
    factor = phase * basis.norms[idx_a] / basis.norms[idx_b]
    return (idx_b, factor)
end

"""
    expand_state(basis::SymmetryReducedBasis{N}, coeffs::AbstractVector)

Expand a symmetry-reduced state vector back into the parent `SpinBasis`.
"""
function expand_state(basis::SymmetryReducedBasis{N}, coeffs::AbstractVector{<:Number}) where {N}
    @assert length(coeffs) == basis.dim "expected $(basis.dim) coefficients, got $(length(coeffs))"

    expanded = zeros(ComplexF64, basis.parent_basis.dim)
    for (parent_idx, state) in enumerate(basis.parent_basis.states)
        rep_idx = representative_index(basis, state)
        rep_idx == 0 && continue
        expanded[parent_idx] = coeffs[rep_idx] * get(basis.state_phases, state, 0.0 + 0.0im) / basis.norms[rep_idx]
    end
    return expanded
end

Base.getindex(basis::SymmetryReducedBasis, i::Int) = basis.representatives[i]
Base.length(basis::SymmetryReducedBasis) = basis.dim
