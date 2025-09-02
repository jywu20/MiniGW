using HDF5
using Test
using LinearAlgebra

const TOL_ZERO = 1e-6

exciton_fid = h5open("/pscratch/sd/j/jywu/WTe2-xy-relaxed/4-absorption/eigenvectors.h5")

@testset "Reciprocal metric" begin
    blat = exciton_fid["mf_header/crystal/blat"] |> read # length unit of b vectors 
    bvec = exciton_fid["mf_header/crystal/bvec"] |> read # Reciprocal vectors
    bvec *= blat

    bdot = exciton_fid["mf_header/crystal/bdot"] |> read # metric of the reciprocal space

    @test norm(bdot - bvec' * bvec) < TOL_ZERO
end

@testset "Dipole" begin
    
end