using LinearAlgebra
using Test
include("../src/MiniGW.jl")
using .MiniGW

const TOL_LARGE = 1e-6
const TOL_SMALL = 1e-10

@testset "Pseudobands" begin
    wfn_full = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
    wfn_pseudobands = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFNo.h5")
    
    # The list comes from /pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/pseudobands.out
    # The full list is too long; 
    # it has been truncated to only keep the first several blocks; 
    # also note that the convention here is the 0-based one.
    #pseudobands_blocks_py_indices_list = [(896, 921), (922, 947), (948, 975), (976, 1003), (1004, 1031), (1032, 1061), (1062, 1091), (1092, 1123), (1124, 1155), (1156, 1189), (1190, 1223), (1224, 1259), (1260, 1295), (1296, 1333), (1334, 1371), (1372, 1411), (1412, 1453), (1454, 1499), (1500, 1545), (1546, 1591), (1592, 1637), (1638, 1685), (1686, 1735), (1736, 1787), (1788, 1839), (1840, 1893), (1894, 1949), (1950, 2007), (2008, 2067), (2068, 2129), (2130, 2193), (2194, 2257), (2258, 2323), (2324, 2393), (2394, 2465), (2466, 2539), (2540, 2615), (2616, 2693), (2694, 2773), (2774, 2857), (2858, 2943), (2944, 3031), (3032, 3121), (3122, 3211), (3212, 3305), (3306, 3403), (3404, 3503), (3504, 3607), (3608, 3713), (3714, 3823), (3824, 3939)]
    pseudobands_blocks_py_indices_list = [(896, 921), (922, 947), (948, 975), (976, 1003), (1004, 1031)]
    last_protected_band_idx = pseudobands_blocks_py_indices_list[1][1] 
    for (block_idx, block_boundary_py) in enumerate(pseudobands_blocks_py_indices_list)
        println("Testing pseudoband #$block_idx")
        pseudoband_idx = last_protected_band_idx + block_idx

        block = block_boundary_py[1] + 1 : block_boundary_py[2] + 1

        ψ_in_block = read_wavefunctions(wfn_full, block, 1 : wfn_full.nrk)
        ψ_pseudo = read_wavefunctions(wfn_pseudobands, pseudoband_idx, 1 : wfn_full.nrk)

        @test sum(ψ_in_block, dims = 3) == ψ_pseudo 
    end
end