using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_range = 1:120
n′ = 4000 
n′′ = 4000 
G_idx = 1000 
k_idx = 12 
q_idx = 37

terms = ComplexF64[]
let progress = Progress(length(n_range), barglyphs=BarGlyphs("[=> ]"), color=:white)
    for n in n_range
        push!(terms, 
            (transition_matrix_irreducible_1BZ(wfn, n, n′,  k_idx, q_idx, G_idx)' * 
            transition_matrix_irreducible_1BZ(wfn, n, n′′, k_idx, q_idx, G_idx))[1])
        next!(progress)
    end
end

let p = plot(dpi=500)
    plot!(p, n_range, real.(terms), framestyle = :box, label = "real")
    plot!(p, n_range, imag.(terms), framestyle = :box, label = "imag")
    
    savefig(p, "nc-n1-$n′-n2-$n′′-nv-$(minimum(n_range))-$(maximum(n_range))-k_idx-$k_idx-q_idx-$q_idx-G_idx-$G_idx.png")
    p
end