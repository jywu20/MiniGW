using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_range = 110:120
n′_range = 2000:2500
G_idx = 200 
k_idx = 12 
q_idx = 37

@time M_nn′ = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, k_idx, q_idx, G_idx)

let p = heatmap(
    n′_range, n′_range, 
    norm.(M_nn′' * M_nn′), 
    aspect_ratio = :equal,
    dpi = 500)
    
    plot_range = (0.5 + first(n′_range), 0.5 + last(n′_range))
    xlims!(p, plot_range)
    ylims!(p, plot_range)
    savefig(p, 
        "nc_range-$(first(n′_range))-$(last(n′_range))-k_idx-$k_idx-q_idx-$q_idx-G_idx-$G_idx.png")
    p
end

begin
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
end