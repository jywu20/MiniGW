using LinearAlgebra
using Test
include("../src/MiniGW.jl")
using .MiniGW
using ProgressMeter

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

"""
Here we only compare the accelerated version 
with the original version; 
if the latter is wrong, 
    then the former is wrong as well, 
    but in a way consistent with the latter.
"""
@testset "BerkeleyGW: accelerated transition matrix: consistency with the definition" begin
    wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    n_range = 117:120
    n′_range = 3999:4000
    k_idx = 3
    q_idx = 4
    G_idx = 200

    M_nn′_kqG = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, 
        k_idx, q_idx, G_idx)
    progress = Progress(length(n_range) * length(n′_range), barglyphs=BarGlyphs("[=> ]"))
    for (new_n, n) in enumerate(n_range)
        for (new_n′, n′) in enumerate(n′_range)
            @test transition_matrix_irreducible_1BZ_def(wfn, n, n′, k_idx, q_idx, G_idx) == M_nn′_kqG[new_n, new_n′]
            next!(progress)
        end
    end
end