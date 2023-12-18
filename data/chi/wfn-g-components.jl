using Plots
using ProgressMeter

include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_idx_rng = 2000:2050
kpt_idx = 68 
σ_idx = 2

G_components_amps = map(n_idx_rng) do n_idx
    ψ = read_wavefunction(wfn, n_idx, kpt_idx)[:, σ_idx]
    G_components = top_n_positions(ψ, 20)
    G_components, abs.(ψ[G_components])
end

n_new_idx = 1
wfn.gvecs_list[kpt_idx][G_components_amps[n_new_idx][1]]
G_components_amps[n_new_idx][2]