using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_range = 60:120
n′_range = 2000:2100
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
        "nc_range-$(first(n′_range))-$(last(n′_range))-nv_range-$(first(n_range))-$(last(n_range))-k_idx-$k_idx-q_idx-$q_idx-G_idx-$G_idx.png")
    p
end