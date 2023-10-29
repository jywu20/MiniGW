using LinearAlgebra
using Test
include("../src/MiniGW.jl")
using .MiniGW

const TOL_LARGE = 1e-6
const TOL_SMALL = 1e-10

@testset "BerkeleyGW: wave function normalization" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    @test abs(norm(read_wavefunction(wfn, 4000, 100)) - 1.0) < TOL_SMALL     
    @test abs(norm(read_wavefunction(wfn, 120,  12))  - 1.0) < TOL_SMALL     
end

@testset "BerkeleyGW: wave function G vectors" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    @test wfn.ngkmax == maximum(wfn.ngk) 
    @test wfn.ngk == length.(wfn.gk_ranges)
end

@testset "BerkeleyGW: wave function individual component" begin
    # TODO
end

@testset "BerkeleyGW: 1BZ check" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    @test equivalent_grids(wfn.full_1BZ, hcat([[x, y, 0.0] for x in -0.45:0.05:0.5 for y in -0.45:0.05:0.5]...))
end

@testset "BerkeleyGW: lattice vector arithmetics" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    @test equivalent_vectors(
        wfn.irreducible_1BZ[:, MiniGW.find_k_plus_q_irreducible_1BZ(wfn, 1, 2)],
        wfn.irreducible_1BZ[:, 1] + wfn.irreducible_1BZ[:, 2]
    ) 

    let k_idx = 10, 
        k_plus_q_idx = 12,
        G_idx = 20,
        G′_idx = 100

        @test equivalent_vectors(
            G_vec(wfn, k_idx, G_idx) + G_vec(wfn, k_idx, G′_idx),
            G_vec(wfn, k_plus_q_idx, MiniGW.find_G_plus_G′(wfn, k_idx, k_plus_q_idx, G_idx, G′_idx))
        ) 
    end
end

# Finding whether G_vec and G_vec_def are consistent 
@testset "BerkeleyGW: finding G vectors" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    k_idx = 120 
    for G_idx in 1 : wfn.ngk[k_idx]
        @test G_vec(wfn, k_idx, G_idx) == MiniGW.G_vec_def(wfn, k_idx, G_idx)
    end
    
    println("Time cost of G_vec:")
    @time for G_idx in 1 : wfn.ngk[k_idx]
        G_vec(wfn, k_idx, G_idx)
    end

    println("Time cost of G_vec_def:")
    @time for G_idx in 1 : wfn.ngk[k_idx]
        MiniGW.G_vec_def(wfn, k_idx, G_idx)
    end
end

# Finding whether find_G_plus_G′ and find_G_plus_G′_def
# are consistent 
@testset "BerkeleyGW: finding G plus G′" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    k_idx = 120 
    k_plus_q_idx = 87
    G′_idx = 145
    G_idx = 2000
    @test MiniGW.find_G_plus_G′(wfn, k_idx, k_plus_q_idx, G_idx, G′_idx) == 
            MiniGW.find_G_plus_G′_def(wfn, k_idx, k_plus_q_idx, G_idx, G′_idx)
end