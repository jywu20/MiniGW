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
G_idx_range = 200:250
k_idx = 1 
q_idx = 1

MM_sum = let res = zeros(ComplexF64, length(n′′_range), length(n′′_range)),
    progress = Progress(length(G_idx_range), barglyphs=BarGlyphs("[=> ]"), color=:white)
    for G_idx in G_idx_range
        M_n′′n_G  = transition_matrix_irreducible_1BZ(wfn, n′′_range, n,  k_idx, q_idx, G_idx)
        M_n′′n′_G = transition_matrix_irreducible_1BZ(wfn, n′′_range, n′, k_idx, q_idx, G_idx)
        res += M_n′′n′_G * M_n′′n_G'
        next!(progress)
    end
    res
end


let p = heatmap(
    n′′_range, n′′_range, 
    norm.(MM_sum), 
    aspect_ratio = :equal,
    dpi = 500)
    
    plot_range = (0.5 + first(n′′_range), 0.5 + last(n′′_range))
    xlims!(p, plot_range)
    ylims!(p, plot_range)
    savefig(p, 
        "nc_range-$(first(n′′_range))-$(last(n′′_range))-n_idx-$n-n_prime_idx_$n′-k_idx-$k_idx-q_idx-$q_idx-G_range-$(minimum(G_idx_range))-$(maximum(G_idx_range)).png")
    p
end

begin
    println("Without pseudobands : $(tr(MM_sum))")
    println("With pseudobands    : $(tr(MM_sum * ones(length(n′′_range), length(n′′_range))))")
end