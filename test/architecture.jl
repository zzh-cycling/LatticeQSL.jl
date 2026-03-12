using Test
using LatticeQSL
using LinearAlgebra

@testset "finite lattice construction" begin
    fl = FiniteLattice(HoneycombLattice(), (2, 3))
    @test fl.nsites == 12
    @test length(fl.bonds) == 18
    @test Set(all_bond_labels(fl)) == Set([:x, :y, :z])
    @test length(bonds_by_label(fl, :x)) == 6
    @test length(bonds_by_label(fl, :y)) == 6
    @test length(bonds_by_label(fl, :z)) == 6
end

@testset "spin basis construction" begin
    basis_full = SpinBasis(4)
    @test length(basis_full) == 16

    basis_sz0 = SpinBasis(4, 2)
    @test length(basis_sz0) == 6
    @test all(num_up(state) == 2 for state in basis_sz0)
end

@testset "heisenberg chain spectrum" begin
    fl = FiniteLattice(ChainLattice(), (2,); pbc = (false,))
    basis = SpinBasis(2, 1)
    model = HeisenbergModel(J = 1.0)
    H = build_hamiltonian(fl, basis, model)
    vals, _ = exact_diagonalization(H)
    @test sort(real.(vals)) ≈ [-0.75, 0.25]
end

@testset "translation reduced basis builds" begin
    fl = FiniteLattice(ChainLattice(), (4,))
    basis = SpinBasis(4, 2)
    tsym = TranslationSymmetry(fl, 1)
    rbasis = SymmetryReducedBasis(basis, tsym, fl)
    @test rbasis.dim > 0
    @test rbasis.dim <= basis.dim

    W = hcat([
        expand_state(rbasis, ComplexF64[ifelse(j == i, 1.0, 0.0) for j in 1:rbasis.dim])
        for i in 1:rbasis.dim
    ]...)
    @test W' * W ≈ Matrix{ComplexF64}(I, rbasis.dim, rbasis.dim)
end

@testset "reflection reduced basis builds" begin
    fl = FiniteLattice(ChainLattice(), (4,))
    basis = SpinBasis(4, 2)
    tsym = TranslationSymmetry(fl, 1)
    rbasis = SymmetryReducedBasis(basis, tsym, fl, 1)
    @test rbasis.dim > 0
    @test rbasis.dim <= basis.dim
end

@testset "entanglement entropy from constrained basis" begin
    basis = SpinBasis(2, 1)
    ψ = ComplexF64[inv(sqrt(2)), -inv(sqrt(2))]
    ρ = reduced_density_matrix(basis, [1], ψ)
    @test ρ ≈ ComplexF64[0.5 0.0; 0.0 0.5]
    @test entanglement_entropy(ρ) ≈ log(2)
end
