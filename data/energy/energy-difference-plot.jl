using LinearAlgebra
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

let n_cs = [1000, 2000, 3000, 4000], n_v = 120
    p = plot(dpi=500) 
    for n_c in n_cs
        plot!(p, 1 : wfn.nrk, wfn.el[n_c, :] - wfn.el[n_v, :], 
        label=n_c,
        framestyle = :box,
        ylims = (0, 1.1 * maximum(wfn.el))) 
    end
    savefig(p, "energies-different-blocks.png")
    p
end

let n_cs = 3825:3940, n_v = 120
    p = plot(dpi=500) 
    for n_c in n_cs
        plot!(p, 1 : wfn.nrk, wfn.el[n_c, :] - wfn.el[n_v, :], 
        legend = false,
        framestyle = :box,
        ylims = (0, 1.1 * maximum(wfn.el))) 
    end
    savefig(p, "energies-same-blocks-$(maximum(n_cs)).png")
    p
end

let n_cs = [125, 126, collect(1000:1200)...], n_v = 120
    p = plot(dpi=500, legend = false) 
    for n_c in n_cs
        plot!(p, 1 : wfn.nrk, 1 ./ (wfn.el[n_c, :] - wfn.el[n_v, :]), 
        framestyle = :box,
    ) 
    end
    p
end
