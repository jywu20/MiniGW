using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_range = 1:120
n′_range = 3939:4000
G_idx = 3000
k_idx = 1 
q_idx = 1

@time M_nn′ = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, k_idx, q_idx, G_idx)

display(M_nn′' * M_nn′)
let p = heatmap(norm.(M_nn′' * M_nn′), aspect_ratio = :equal)
    xlims!(p, (0.5, length(n′_range) + 0.5))
    ylims!(p, (0.5, length(n′_range) + 0.5))
    p
end

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