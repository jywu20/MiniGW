using LinearAlgebra
using Test
using ProgressMeter
include("../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_range = 117:120
n′_range = 1999:2002
G_idx = 810 
k_idx = 2
q_idx = 3

M_nn′ = zeros(ComplexF64, length(n_range), length(n′_range))

let progress = Progress(length(M_nn′), barglyphs=BarGlyphs("[=> ]"))
    for (n_idx, n) in enumerate(n_range)
        for (n′_idx, n′) in enumerate(n′_range) 
            M_nn′[n_idx, n′_idx] = 
                transition_matrix_irreducible_1BZ_def(wfn, n, n′, k_idx, q_idx, G_idx)
            next!(progress)
        end
    end
end

display(M_nn′)

let χ = 0.0
    for n in eachindex(n_range)
        for n′ in eachindex(n′_range)
            χ += M_nn′[n, n′]' * M_nn′[n, n′]
        end
    end

    println("With no pseudobands : $χ")
end

let χ_pb = 0.0 
    for n in eachindex(n_range)
        for n′_1 in eachindex(n′_range)  
            for n′_2 in eachindex(n′_range)
                χ_pb += M_nn′[n, n′_1]' * M_nn′[n, n′_2]
            end
        end
    end
    
    println("With pseudobands    : $χ_pb")
end