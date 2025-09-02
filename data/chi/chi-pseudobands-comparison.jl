using LinearAlgebra
using Statistics
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

# What happens when valance bands are pseudolized?
let k_idx_range = 1:121,
    q_idx = 1,
    n_1 = 118,
    n_2 = 118,
    n′ = 140,
    G_idx = 100,
    G′_idx = 102
    
    if n_1 > n_2
        temp = n_1
        n_1 = n_2
        n_2 = n_1
    end
    
    χ_terms = ComplexF64[] 
    progress = Progress(length(k_idx_range), barglyphs=BarGlyphs("[=> ]"), color=:white) 
    for k_idx in k_idx_range
        inverse_ΔE = 1 / (mean(wfn.el[n_1:n_2, k_idx_range]) - wfn.el[n′, k_idx])

        M_n1n′_kqG  = transition_matrix_irreducible_1BZ(wfn, n_1, n′, k_idx, q_idx, G_idx)[1]
        M_n2n′_kqG′ = transition_matrix_irreducible_1BZ(wfn, n_2, n′, k_idx, q_idx, G′_idx)[1]


        push!(χ_terms, M_n1n′_kqG * M_n2n′_kqG′' * inverse_ΔE)  
        next!(progress)
    end
    
    println(sum(χ_terms))
    p = plot()
    plot!(p, real.(χ_terms)) 
    plot!(p, imag.(χ_terms))
    p
end

# What happens when valance bands are pseudolized?
let k_idx_range = 1:121,
    q_idx = 1,
    n_1 = 118,
    n_2 = 118,
    n′ = 121,
    G_idx = 100,
    G′_idx = 102
    
    if n_1 > n_2
        temp = n_1
        n_1 = n_2
        n_2 = n_1
    end
    
    χ_terms = ComplexF64[] 
    progress = Progress(length(k_idx_range), barglyphs=BarGlyphs("[=> ]"), color=:white) 
    for k_idx in k_idx_range
        inverse_ΔE = 1 / (mean(wfn.el[n_1:n_2, k_idx]) - wfn.el[n′, k_idx])

        M_n1n′_kqG  = transition_matrix_irreducible_1BZ(wfn, n_1, n′, k_idx, q_idx, G_idx)[1]
        M_n2n′_kqG′ = transition_matrix_irreducible_1BZ(wfn, n_2, n′, k_idx, q_idx, G′_idx)[1]


        push!(χ_terms, M_n1n′_kqG * M_n2n′_kqG′' * inverse_ΔE)  
        next!(progress)
    end
    
    println(sum(χ_terms))
    p = plot()
    plot!(p, real.(χ_terms)) 
    plot!(p, imag.(χ_terms))
    p
end