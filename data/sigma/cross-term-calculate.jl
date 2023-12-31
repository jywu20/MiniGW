using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n = 120 
n′ = 120
n′′_range = 1800:2000
G_idx = 200
k_idx = 1 
q_idx = 1

@time M_n′′n_G  = transition_matrix_irreducible_1BZ(wfn, n′′_range, n,  k_idx, q_idx, G_idx)
@time M_n′′n′_G = transition_matrix_irreducible_1BZ(wfn, n′′_range, n′, k_idx, q_idx, G_idx)

let p = heatmap(
    n′′_range, n′′_range, 
    norm.(M_n′′n′_G * M_n′′n_G'), 
    aspect_ratio = :equal,
    dpi = 500)
    
    plot_range = (0.5 + first(n′′_range), 0.5 + last(n′′_range))
    xlims!(p, plot_range)
    ylims!(p, plot_range)
    savefig(p, 
        "nc_range-$(first(n′′_range))-$(last(n′′_range))-n_idx-$n-n_prime_idx_$n′-k_idx-$k_idx-q_idx-$q_idx-G_idx-$G_idx.png")
    p
end
