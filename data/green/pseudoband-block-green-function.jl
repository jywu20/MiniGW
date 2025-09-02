using LinearAlgebra
using Plots
using Colors
include("../../src/MiniGW.jl")
using .MiniGW
using Statistics

wfn_full = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

#pseudobands_blocks_py_indices_list = [ (3608, 3713), (3714, 3823), (3824, 3939)]
pseudobands_blocks_py_indices_list = [ (3608, 3650), (3651, 3713), (3714, 3770), (3771, 3823), (3824, 3880), (3881, 3939)]

bottom_idx = pseudobands_blocks_py_indices_list[1][1] + 1
top_idx    = pseudobands_blocks_py_indices_list[end][2] + 1
kpt_idx = 1

N_ξ = 10

function projector(bottom_idx, top_idx, n)
    total_subspace_size = top_idx - bottom_idx + 1
    P = zeros(ComplexF64, total_subspace_size, total_subspace_size)
    P[n - bottom_idx + 1, n - bottom_idx + 1] = 1.0 
    P
end

function band_basis(bottom_idx, top_idx, n)
    total_subspace_size = top_idx - bottom_idx + 1
    ψ = zeros(ComplexF64, total_subspace_size)
    ψ[n - bottom_idx + 1] = 1.0
    ψ
end

G_correct = sum(bottom_idx : top_idx) do n
    projector(bottom_idx, top_idx, n) / read_energies(wfn_full, n, kpt_idx)
end

G_pseudo = sum(pseudobands_blocks_py_indices_list) do block
    block_subspace = block[1] + 1 : block[2] + 1
    sum(1 : N_ξ) do ξ_idx
        ξ = sum(block_subspace) do band_idx
            exp(2π * im * rand()) * band_basis(bottom_idx, top_idx, band_idx)
        end / sqrt(N_ξ)
        
        E = mean(read_energies(wfn_full, block_subspace, :))
        ξ * ξ' / E
    end
end

heatmap(bottom_idx : top_idx, bottom_idx : top_idx, 
    abs.((G_correct - G_pseudo) / tr(G_correct) * (top_idx - bottom_idx + 1)),
    c = :bilbao50, 
    aspect_ratio = :equal, 
    xlims = (bottom_idx - 0.5, top_idx + 0.5), 
    ylims = (bottom_idx - 0.5, top_idx + 0.5), 
    frame_style = :box)

# The following code is used to check whether the non-diagonal elements in G_pseudo 
# indeed have reasonable values, 
# when the 1/E factor is not included; 
# it seems the answer is "yes".
let N_ξ = 40, arr = zeros(ComplexF64, N_ξ) 
    for i = 1 : 40
        α1 = exp(2π * im * rand())
        α2 = exp(2π * im * rand())
        arr[i] = α1' * α2
    end
    
    plot(1 : N_ξ, map(1 : N_ξ) do N_ξ′
        abs(sum(arr[1 : N_ξ′])) / N_ξ′
    end)
end

# Anyway, what is known for sure is that the parameters used in pseudobands.py 
# can't lead to fully converged DFT Green function.