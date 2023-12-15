using LinearAlgebra
using Plots
using Colors
include("../../src/MiniGW.jl")
using .MiniGW
using Statistics

wfn_full = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
wfn_pseudobands = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFNo.h5")

pseudobands_blocks_py_indices_list = [(896, 921), (922, 947), (948, 975), (976, 1003), (1004, 1031), (1032, 1061), (1062, 1091), (1092, 1123), (1124, 1155), (1156, 1189), (1190, 1223), (1224, 1259), (1260, 1295), (1296, 1333), (1334, 1371), (1372, 1411), (1412, 1453), (1454, 1499), (1500, 1545), (1546, 1591), (1592, 1637), (1638, 1685), (1686, 1735), (1736, 1787), (1788, 1839), (1840, 1893), (1894, 1949), (1950, 2007), (2008, 2067), (2068, 2129), (2130, 2193), (2194, 2257), (2258, 2323), (2324, 2393), (2394, 2465), (2466, 2539), (2540, 2615), (2616, 2693), (2694, 2773), (2774, 2857), (2858, 2943), (2944, 3031), (3032, 3121), (3122, 3211), (3212, 3305), (3306, 3403), (3404, 3503), (3504, 3607), (3608, 3713), (3714, 3823), (3824, 3939)]

p = plot(legend = false)

let first_idx = pseudobands_blocks_py_indices_list[1][1] + 1, 
    last_idx  = pseudobands_blocks_py_indices_list[end][2] + 1
    
    high_bands = first_idx : last_idx
    scatter!(p, high_bands, 1 ./ read_energies(wfn_full, high_bands, 1), 
    markerstrokewidth=0)
end

pseudobands_blocks_E_inv = zeros(length(pseudobands_blocks_py_indices_list))
for (idx_block, block) in enumerate(pseudobands_blocks_py_indices_list) 
    block_indices = block[1] + 1 : block[2] + 1
    el = read_energies(wfn_full, block_indices, :) 
    mean_E_inv = mean(1 ./ el)
    pseudobands_blocks_E_inv[idx_block] = mean_E_inv
    plot!(p, block_indices, mean_E_inv * ones(length(block_indices)), c = :salmon1)
end

p

# It can be observed that for each energy region, 
# we have several pseudoband blocks, 
# which likely is enough to guarantee the stochastic convergence