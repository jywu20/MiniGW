using Plots
using ProgressMeter

include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")


function MM′_GG′(G_1_idx, G_2_idx, nv_range, k_idx, σ_idx)
    map(nv_range) do nv
        read_wavefunction(wfn, nv, k_idx)[G_1_idx, σ_idx]' * 
        read_wavefunction(wfn, nv, k_idx)[G_2_idx, σ_idx]
    end
end

# Verify that ∑_n^occ c_n(G_1) c_n(G_2) vanishes.
let G_1_idx = 1940,
    G_2_idx = 1960,
    nv_range = 1 : 120

    MM′ = MM′_GG′(G_1_idx, G_2_idx, nv_range, 56, 1)
    
    println(sum(MM′))
    plot(real.(MM′)) 
end

# Verify that ∑_n^occ c_n(G_1) c_n(G_2) vanishes, 
# for a wide range of G vectors
let G_list = 1940:1960, 
    nv_range = 1 : 120, 
    k_idx = 2, 
    σ_idx = 2
    progress = Progress(length(G_list)^2)
    MM′_GG′_list = map(Iterators.product(G_list, G_list)) do (G_1_idx, G_2_idx)
        next!(progress)
        abs(sum(MM′_GG′(G_1_idx, G_2_idx, nv_range, k_idx, σ_idx))) 
    end

    heatmap(G_list, G_list, MM′_GG′_list', 
        aspect_ratio = :equal, 
        xlims = (minimum(G_list) + 0.5, maximum(G_list) - 0.5),
        ylims = (minimum(G_list) + 0.5, maximum(G_list) - 0.5),
    )
end

# Verify that G vectors connected by symmetry lead to same MM'
let G_idx_list = [1944, 1947], 
    σ_idx = 1
    sum(1 : 120) do nv
        abs.(read_wavefunction(wfn, nv, 1)[G_idx_list, σ_idx]).^2
    end    
end

# Pseudobands for valence bands
let nv_range = 1 : 120,
    G_1_idx = 100,
    G_2_idx = 100,
    kpt_idx = 12, 
    σ_idx = 1 
    
    progress = Progress(length(nv_range)^2)
    sum(Iterators.product(nv_range, nv_range)) do (n_1, n_2) 
        next!(progress)
        read_wavefunction(wfn, n_1, kpt_idx)[G_1_idx, σ_idx]' * 
        read_wavefunction(wfn, n_2, kpt_idx)[G_2_idx, σ_idx]
    end |> println
    
    println(
        sum(nv_range) do n_1
            read_wavefunction(wfn, n_1, kpt_idx)[G_1_idx, σ_idx]' * 
            read_wavefunction(wfn, n_1, kpt_idx)[G_2_idx, σ_idx]
        end
    )
end