using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_range = 1:120
n′_range = 1951:2000
G_idx = 810 
G′_idx = 1000 
k_idx = 2
q_idx = 3

@time M_nn′_G  = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, k_idx, q_idx, G_idx)
@time M_nn′_G′ = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, k_idx, q_idx, G′_idx)


let p = heatmap(
    n′_range, n′_range, 
    norm.(M_nn′_G' * M_nn′_G′), 
    aspect_ratio = :equal,
    dpi = 500)
    
    plot_range = (0.5 + first(n′_range), 0.5 + last(n′_range))
    xlims!(p, plot_range)
    ylims!(p, plot_range)
    savefig(p, 
        "nc_range-$(first(n′_range))-$(last(n′_range))-k_idx-$k_idx-q_idx-$q_idx-G1_idx-$G_idx-G2_idx-$G′_idx.png")
    p
end

begin
    let χ = 0.0
        for n in eachindex(n_range)
            for n′ in eachindex(n′_range)
                χ += M_nn′_G[n, n′]' * M_nn′_G′[n, n′]
            end
        end
    
        println("With no pseudobands : $χ")
    end
    
    let χ_pb = 0.0 
        for n in eachindex(n_range)
            for n′_1 in eachindex(n′_range)  
                for n′_2 in eachindex(n′_range)
                    χ_pb += M_nn′_G[n, n′_1]' * M_nn′_G′[n, n′_2]
                end
            end
        end
        
        println("With pseudobands    : $χ_pb")
    end
end

let p = heatmap(
    n′_range, n′_range, 
    real.(M_nn′_G' * M_nn′_G′), 
    aspect_ratio = :equal,
    dpi = 500)
    
    plot_range = (0.5 + first(n′_range), 0.5 + last(n′_range))
    xlims!(p, plot_range)
    ylims!(p, plot_range)
    savefig(p, 
        "nc_range-$(first(n′_range))-$(last(n′_range))-k_idx-$k_idx-q_idx-$q_idx-G1_idx-$G_idx-G2_idx-$G′_idx-real-part.png")
    p
end