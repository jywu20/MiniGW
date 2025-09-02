using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n = 120 
n′ = 120
n′′_range = 1800:2000
G_idx_range = 200:250
k_idx = 120 
q_idx = 1

MM_sum = let res = zeros(ComplexF64, length(n′′_range), length(n′′_range)),
    progress = Progress(length(G_idx_range), barglyphs=BarGlyphs("[=> ]"), color=:white)

    M_n′′n_G  = Vector{Matrix{ComplexF64}}(undef, length(G_idx_range)) 
    M_n′′n′_G = Vector{Matrix{ComplexF64}}(undef, length(G_idx_range)) 
    for (G_idx_new, G_idx) in enumerate(G_idx_range) 
        M_n′′n_G[G_idx_new]  = transition_matrix_irreducible_1BZ(wfn, n′′_range, n,  k_idx, q_idx, G_idx)
        M_n′′n′_G[G_idx_new] = transition_matrix_irreducible_1BZ(wfn, n′′_range, n′, k_idx, q_idx, G_idx)
        next!(progress)
    end
    
    for G_idx_new in eachindex(G_idx_range)
        for G′_idx_new in eachindex(G_idx_range)
            res += M_n′′n_G[G_idx_new] * M_n′′n′_G[G′_idx_new]'
        end
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
        "nc_range-$(first(n′′_range))-$(last(n′′_range))-n_idx-$n-n_prime_idx_$n′-k_idx-$k_idx-q_idx-$q_idx-G_range-$(minimum(G_idx_range))-$(maximum(G_idx_range))-double-G.png")
    p
end

begin
    println("Without pseudobands : $(tr(MM_sum))")
    println("With pseudobands    : $(tr(MM_sum * ones(length(n′′_range), length(n′′_range))))")
end