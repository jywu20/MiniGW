using LinearAlgebra
using Test
include("../src/MiniGW.jl")
using .MiniGW

const TOL_LARGE = 1e-6
const TOL_SMALL = 1e-10

@testset "Utilities: grid-related" begin
    grid = [
        1 2 3;
        4 5 6
    ]
    @test find_in_grid(grid, [2, 5]) == 2
    @test find_in_grid(grid, [2, 6]) == -1
end

@testset "Utilities: periodic grid-related" begin
    @test abs(wigner_seitz(1.4) - 0.4) < TOL_SMALL
    @test abs(wigner_seitz(-0.6) - 0.4) < TOL_SMALL
    @test abs(wigner_seitz(0.49) - 0.49) < TOL_SMALL
    @test abs(wigner_seitz(2.49) - 0.49) < TOL_SMALL
    @test abs(wigner_seitz(0.5) - 0.5) < TOL_SMALL
    @test abs(wigner_seitz(-0.5) - 0.5) < TOL_SMALL
    
    grid = [
        0.5  0.6  0.7  0.4; 
        0.0  0.1  0.2  0.3;
        0.1 -0.1 -0.3  1.3
    ]
    @test find_in_periodic_grid(grid, [0.4, 0.3, 1.3]) == 4
    @test find_in_periodic_grid(grid, [0.4, 0.3, 1.2]) == -1
    @test find_in_periodic_grid(grid, [1.5, -1.0, -10.9]) == 1
end

@testset "Utilities: symmetry expansion of 1BZ" begin
    irreducible_1BZ = [
        0.0 0.1  
        0.0 0.0
    ]
    ops = [
        [1 0; 0 1], 
        [-1 0; 0 -1]
    ]
    @test expand_sym_reduced_grid(irreducible_1BZ, ops) == [
        0.0 0.1 -0.1
        0.0 0.0  0.0
    ]

    irreducible_1BZ = hcat([0.1 * [i, j, 0.0] for i in 0 : 5 for j in 0 : 5]...)
    ops = [
        [
            1 0 0; 
            0 1 0;
            0 0 1
        ],
        [
            -1  0  0;
             0 -1  0;
             0  0  1
        ],
        [
             1  0  0;
             0 -1  0;
             0  0  1
        ],
        [
            -1  0  0;
             0  1  0;
             0  0  1
        ],
    ]
    full_1BZ = hcat(
        [0.1 * [i, j, 0.0] for i in -4 : 5 for j in -4 : 5]...)
    @test equivalent_grids(expand_sym_reduced_grid(irreducible_1BZ, ops), full_1BZ) 
end
