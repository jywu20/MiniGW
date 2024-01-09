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

# Bohr radius, in Å 
a0 = 0.5291772
# Crystal constants, again in Å
a = 6.314769
b = 3.492485
c = 20.16

kx, ky, kz = wfn.irreducible_1BZ[:, kpt_idx]

function Ecut(gvec)
    Gx, Gy, Gz = gvec
    (2π * a0 / a)^2 * (kx + Gx)^2 + 
    (2π * a0 / b)^2 * (ky + Gy)^2 + 
    (2π * a0 / c)^2 * (kz + Gz)^2
end

n_new_idx = 1
high_amp_G_vec_idx_rng = G_components_amps[n_new_idx][1] 
high_amp_G_vecs = wfn.gvecs_list[kpt_idx][high_amp_G_vec_idx_rng]
G_amps = G_components_amps[n_new_idx][2]

E_cut_list = map(Ecut, wfn.gvecs_list[kpt_idx][G_components_amps[n_new_idx][1]])

#G = [3, 2, 5]
G = [0, 0, 0]
G_idx = find_in_grid(wfn.gvecs[kpt_idx], G)

map(high_amp_G_vec_idx_rng) do G′_idx
    G_plus_G′_idx = MiniGW.find_G_plus_G′(wfn, kpt_idx, kpt_idx, G_idx, G′_idx)
    abs(read_wavefunction(wfn, 120, kpt_idx)[G_plus_G′_idx, σ_idx])^2
end
# Large dispersion can still be observed hmm.