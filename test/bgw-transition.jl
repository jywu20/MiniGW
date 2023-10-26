using LinearAlgebra
using Test
include("../src/MiniGW.jl")
using .MiniGW

wfn = BerkeleyGWSpinorWaveFunction("/pscratch/sd/j/jywu/WTe2-xy-relaxed/2.1-wfn-xy/WFN.h5")
@time M_nn′_kqG = transition_matrix_def(wfn, 120, 122, 2, 1, 2)
println(M_nn′_kqG)