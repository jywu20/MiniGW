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
        wfn.irreducible_1BZ[:, MiniGW.find_k_plus_q(wfn, 1, 2)],
        wfn.irreducible_1BZ[:, 1] + wfn.irreducible_1BZ[:, 2]
    ) 

    let k_idx = 10, 
        k_plus_q_idx = 12,
        G_idx = 20,
        G′_idx = 100

        @test equivalent_vectors(
            wfn.gvecs[k_idx][:, G_idx] + wfn.gvecs[k_idx][:, G′_idx],
            wfn.gvecs[k_plus_q_idx][:, MiniGW.find_G_plus_G′(wfn, k_idx, k_plus_q_idx, G_idx, G′_idx)]
        ) 
    end
end