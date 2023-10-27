using LinearAlgebra
using Test
include("../src/MiniGW.jl")
using .MiniGW

@testset "BerkeleyGW: definition of transition matrix: normalization" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

    # When q = 0, G = 0 and n and n′ are the same,  
    # the transition matrix should always be 1.
    @time M_nn′_kqG = transition_matrix_irreducible_1BZ_def(wfn, 120, 120, 2, 1, 1)
    @test isapprox(M_nn′_kqG, 1.0)
    @time M_nn′_kqG = transition_matrix_irreducible_1BZ_def(wfn, 122, 122, 40, 1, 1)
    @test isapprox(M_nn′_kqG, 1.0)
    @time M_nn′_kqG = transition_matrix_irreducible_1BZ_def(wfn, 4000, 4000, 80, 1, 1)
    @test isapprox(M_nn′_kqG, 1.0)
end
