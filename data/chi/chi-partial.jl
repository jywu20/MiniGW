using LinearAlgebra
using Statistics
using Test
using ProgressMeter
using Plots
include("../../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWFN("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")

n_range = 1:120
G_idx = 800 
k_idx = 15 
q_idx = 80

# Read from /pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/pseudobands.out
pseudobands_blocks = [(2194, 2257), (2258, 2323), (2324, 2393)]

#n′_range = 2000:2050
n′_range = (pseudobands_blocks[1][1]:pseudobands_blocks[end][2]) .+ 1

# Not pseudobands
let χ = 0.0
    @time M_nn′_kqG = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, k_idx, q_idx, G_idx)
    for (n′_idx, n′) in enumerate(n′_range) 
        for (n_idx, n) in enumerate(n_range) 
            inverse_ΔE = 1 / transition_energy(wfn, n, n′, k_idx, q_idx)
            χ += M_nn′_kqG[n_idx, n′_idx] * M_nn′_kqG[n_idx, n′_idx]' * inverse_ΔE
        end
    end
    χ
end

# Pseudobands: regard n′_range as a whole pseudoband block 
let χ = 0.0, 
    progress = Progress(length(n′_range), barglyphs=BarGlyphs("[=> ]"), color=:white)
    @time M_nn′_kqG = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, k_idx, q_idx, G_idx)
    E_c_avg = mean(wfn.el[n′_range, :]) 
    for (n′_1_idx, n′) in enumerate(n′_range) 
        for n′_2_idx in 1 : length(n′_range)
            for (n_idx, n) in enumerate(n_range) 
                inverse_ΔE = 1 / (wfn.el[n, k_idx] - E_c_avg) 
                χ += M_nn′_kqG[n_idx, n′_1_idx] * M_nn′_kqG[n_idx, n′_2_idx]' * inverse_ΔE
            end
        end
        next!(progress)
    end
    χ
end

# Pseudobands: use blocks in pseudobands_blocks
let χ = 0.0, 
    progress = Progress(length(pseudobands_blocks), barglyphs=BarGlyphs("[=> ]"), color=:white)
    @time M_nn′_kqG = transition_matrix_irreducible_1BZ(wfn, n_range, n′_range, k_idx, q_idx, G_idx)
    for current_block_original in pseudobands_blocks
        current_block = (current_block_original[1] : current_block_original[2]) .+ 1
        E_c_avg = mean(wfn.el[current_block, :])
        for (n′_1_idx, n′) in enumerate(current_block) 
            for n′_2_idx in 1 : length(current_block)
                for (n_idx, n) in enumerate(n_range) 
                    inverse_ΔE = 1 / (wfn.el[n, k_idx] - E_c_avg) 
                    χ += M_nn′_kqG[n_idx, n′_1_idx] * M_nn′_kqG[n_idx, n′_2_idx]' * inverse_ΔE
                end
            end
        end
        next!(progress)
    end
    χ
end