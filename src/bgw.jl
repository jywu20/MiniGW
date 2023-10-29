export Ry_BGW, 
    BerkeleyGWSpinorWaveFunction,
    G_vec,
    read_wavefunction,
    read_wavefunctions, 
    transition_matrix_irreducible_1BZ_def,
    transition_matrix_irreducible_1BZ

######################################################################
#region IO

const Ry_BGW =  13.6056925

struct BerkeleyGWSpinorWaveFunction <: AbstractGWWaveFunction
    fid::HDF5.File
    
    nrk::Int
    irreducible_1BZ::Matrix{Float64}
    full_1BZ::Matrix{Float64}
    full_to_irreducible::Vector{Int}
    el::Matrix{Float64}

    gk_ranges::Vector{UnitRange} 
    ngk::Vector{Int} 
    ngkmax::Int
    gvecs::Vector{Matrix{Int}}
    gvecs_list::Vector{Vector{SVector{3, Int}}}
end

function BerkeleyGWSpinorWaveFunction(fid::HDF5.File; tol_sym = 1e-6)
    if read(fid["mf_header/kpoints/nspinor"]) != 2
        throw("The HDF5 file is not a spinor wave function.")
    end
    
    if read(fid["mf_header/kpoints/nspin"]) == 2
        throw("The HDF5 file comes from a collinear calculation.")
    end

    el = read(fid["mf_header/kpoints/el"])[:, :, 1] * Ry_BGW

    nrk = fid["mf_header/kpoints/nrk"] |> read
    rk = fid["mf_header/kpoints/rk"] |> read
    
    ntran = fid["mf_header/symmetry/ntran"] |> read
    full_1BZ, full_to_irreducible = 
    expand_sym_reduced_grid(rk, map(1 : ntran) do sym_idx
        fid["mf_header/symmetry/mtrx"][:, :, sym_idx]
    end, tol = tol_sym)
    
    # Finding how many G vectors each k point has 
    ngk = fid["mf_header/kpoints/ngk"] |> read 
    ngkmax = fid["mf_header/kpoints/ngkmax"] |> read
    gk_ending_points = map(eachindex(ngk)) do kpt_idx
        sum(ngk[1 : kpt_idx])
    end
    gk_ending_points = [0, gk_ending_points...]
    gk_ranges = map(eachindex(ngk)) do kpt
        gk_ending_points[kpt] + 1 : gk_ending_points[kpt + 1]
    end
    
    # Finding G vectors for each k point 
    gvecs_original_form = fid["wfns/gvecs"] |> read
    gvecs = Vector{Matrix{Int}}(undef, nrk)
    gvecs_list = Vector{Vector{Vector{Int}}}(undef, nrk)
    for k_idx in 1 : nrk
        gvecs[k_idx] = gvecs_original_form[:, gk_ranges[k_idx]]
        gvecs_list[k_idx] = grid_to_list(gvecs[k_idx])
    end
    
    BerkeleyGWSpinorWaveFunction(fid, 
        nrk, rk, full_1BZ, full_to_irreducible,
        el, 
        gk_ranges, ngk, ngkmax, 
        gvecs, gvecs_list)
end

function BerkeleyGWSpinorWaveFunction(path::AbstractString; allow_write = false)
    if allow_write
        fid = h5open(path, "r+")
    else
        fid = h5open(path)
    end
    BerkeleyGWSpinorWaveFunction(fid)
end

"""
Close the HDF5 file associated to `wfn`; 
this doesn't delete the existing data extracted from `wfn.fid`.
"""
function Base.close(wfn::BerkeleyGWSpinorWaveFunction)
    close(wfn.fid)
end

"""
The convention of the order of array indices is this:
`read_wavefunction(wfn::BerkeleyGWWaveFunction, band_idx, kpt_idx)[G_idx, σ_idx]`.
Note that in the current implementation, 
the `band_idx` can be a range; 
in this case the return value has the following form:
`read_wavefunction(wfn::BerkeleyGWWaveFunction, band_idx, kpt_idx)[G_idx, σ_idx, band_idx]`.
"""
function read_wavefunction(wfn::BerkeleyGWSpinorWaveFunction, 
    band_idx, kpt_idx; extend_G_grid = false)

    fid = wfn.fid
    gk_ranges = wfn.gk_ranges

    Re_c_nk_Gσ = fid["wfns/coeffs"][1, gk_ranges[kpt_idx], :, band_idx]
    Im_c_nk_Gσ = fid["wfns/coeffs"][2, gk_ranges[kpt_idx], :, band_idx]
    
    compelx_wavefunction = Re_c_nk_Gσ + im * Im_c_nk_Gσ
    if ! extend_G_grid
        return compelx_wavefunction
    end 

    # TODO: there seems to be a bug here:
    # the same G vector index may refer to two different G vectors  
    # for two k points, 
    # when we are at the edge of the G grid 
    # This shouldn't be a major headache right now 
    # since the cutoff energy is usually large enough; 
    # still it has to be corrected at a certain stage 
    # Currently the extend_G_grid flag should always be false
    ngkmax = wfn.ngkmax
    extended_complex_wavefunction = zeros(ComplexF64, ngkmax, 2)
    extended_complex_wavefunction[1 : size(compelx_wavefunction)[1], :] = compelx_wavefunction
    extended_complex_wavefunction
end

"""
The convention of the order of array indices is this:
read_wavefunction(wfn::BerkeleyGWWaveFunction, band_idx, kpt_idx)[G_idx, σ_idx, band_idx, kpoint_idx]

Note that `band_idx` and `kpoint_idx` here are indices within `band_range` and `kpt_range`; 
users need to keep a copy of `band_range` and `kpt_range`
to establish the correspondence between the indices in the total wave function file 
and the indices in the block extracted.
"""
function read_wavefunctions(wfn::BerkeleyGWSpinorWaveFunction, 
    band_range, kpt_range)
    
    ngkmax = wfn.ngkmax

    wavefunctions = zeros(ComplexF64, ngkmax, 2, length(band_range), length(kpt_range))
    for (new_band_idx, band_idx) in enumerate(band_range), 
        (new_kpt_idx, kpt_idx) in enumerate(kpt_range) 
        wavefunctions[:, :, new_band_idx, new_kpt_idx] = 
            read_wavefunction(wfn, band_idx, kpt_idx, extend_G_grid = true)
    end
    
    wavefunctions
end

#endregion 
######################################################################

######################################################################
#region k-points and G points 

function G_vec_def(wfn::BerkeleyGWSpinorWaveFunction, k_idx, G_idx)
    wfn.gvecs[k_idx][:, G_idx]
end

function G_vec_list(wfn::BerkeleyGWSpinorWaveFunction, k_idx, G_idx)
    wfn.gvecs_list[k_idx][G_idx]
end

G_vec = G_vec_list

"""
`q_idx` is assumed to be the index in `wfn.irreducible`; 
some optimization may be done here.
"""
function find_k_plus_q_irreducible_1BZ(wfn::BerkeleyGWSpinorWaveFunction, k_idx, q_idx)
    irreducible_1BZ = wfn.irreducible_1BZ
    full_1BZ = wfn.full_1BZ
    k_plus_q =  irreducible_1BZ[:, k_idx] + irreducible_1BZ[:, q_idx]
    k_plus_q_in_full_1BZ = find_in_periodic_grid(full_1BZ, k_plus_q)
    wfn.full_to_irreducible[k_plus_q_in_full_1BZ]
end

"""
`G_idx` is an index in the polarizability matrix; 
usually the GW cutoff energy is chosen to be much smaller 
than the DFT cutoff energy, 
so we assume that in the G grid of `k_idx` and the G grid of `k_plus_q`, 
`G_idx` should refer to the same thing.
`G′_idx` is an index in the G grid of `k_idx`.
Some optimization may be done here.
"""
function find_G_plus_G′_def(wfn::BerkeleyGWSpinorWaveFunction, k_idx, k_plus_q_idx, G_idx, G′_idx)
    G_grid_of_k_plus_q = wfn.gvecs[k_plus_q_idx]
    G  = G_vec(wfn, k_idx, G_idx) 
    G′ = G_vec(wfn, k_idx, G′_idx) 
    find_in_grid(G_grid_of_k_plus_q, G + G′)
end

function find_G_plus_G′(wfn::BerkeleyGWSpinorWaveFunction, k_idx, k_plus_q_idx, G_idx, G′_idx)
    G_grid_of_k_plus_q = wfn.gvecs_list[k_plus_q_idx]
    G  = G_vec(wfn, k_idx, G_idx) 
    G′ = G_vec(wfn, k_idx, G′_idx) 
    find_in_list(G_grid_of_k_plus_q, G + G′)
end

function indices_of_G_plus_G′_def(wfn::BerkeleyGWSpinorWaveFunction, 
    k_idx, q_idx, G_idx)
    k_plus_q_idx = find_k_plus_q_irreducible_1BZ(wfn, k_idx, q_idx)
    
    # Since G′ is confined in the G grid of k, 
    # the number of G + G′, once G is given, 
    # is the same as the size of the G grid of k.
    G_plus_G′_indices = zeros(Int, wfn.ngk[k_idx])

    for G′_idx in 1 : wfn.ngk[k_idx]
        G_plus_G′_indices[G′_idx] = 
            find_G_plus_G′(wfn, k_idx, k_plus_q_idx, G_idx, G′_idx)
    end
    
    G_plus_G′_indices
end

"""
Return an array with the same size of the G grid of `k_idx`,
which is the range of G′; 
the `G′_idx`-th element of the return value 
is the index of G+G′ in the G grid of k+q.
"""
indices_of_G_plus_G′ = indices_of_G_plus_G′_def
# The implementation below seems to also have performance problem:
#function indices_of_G_plus_G′(wfn::BerkeleyGWSpinorWaveFunction, k_idx, k_plus_q_idx, G_idx)
#    G_grid_of_k = wfn.gvecs[k_idx] 
#    G_grid_of_k_plus_q = wfn.gvecs[k_plus_q_idx]
#
#    G = G_grid_of_k[:, G_idx]
#    # By default the return value is -1
#    G_plus_G′_grid = -ones(wfn.ngk[k_plus_q_idx])
#    possible_G_plus_G′_indices = 1 : wfn.ngk[k_plus_q_idx]
#    existing_G_plus_G′_indices = Set{Int}()
#    
#    for G′_idx in 1 : wfn.ngk[k_idx]
#        G′ = G_grid_of_k[:, G′_idx]
#        for possible_G_plus_G′_idx in possible_G_plus_G′_indices
#            if possible_G_plus_G′_idx in existing_G_plus_G′_indices
#                continue
#            end
#
#            if G + G′ == G_grid_of_k_plus_q[:, possible_G_plus_G′_idx]
#                G_plus_G′_grid[G′_idx] = possible_G_plus_G′_idx
#                push!(existing_G_plus_G′_indices, possible_G_plus_G′_idx)
#            end
#        end
#    end
#
#    G_plus_G′_grid 
#end

#endregion
######################################################################

######################################################################
#region Transition matrix 

"""
The most naive implementation of the transition matrix M_nn'(k, q, G).
`n` and `n′` should be scalars.
"""
function transition_matrix_irreducible_1BZ_def(wfn::BerkeleyGWSpinorWaveFunction, 
    n, n′, k_idx, q_idx, G_idx)
    
    k_plus_q_idx = find_k_plus_q_irreducible_1BZ(wfn, k_idx, q_idx)
    
    c_n′k_Gσ       = read_wavefunction(wfn, n′, k_idx)
    c_nk_plus_q_Gσ = read_wavefunction(wfn, n, k_plus_q_idx)

    sum(map(1 : wfn.ngk[k_idx]) do G′_idx
        G_plus_G′_idx = find_G_plus_G′(wfn, k_idx, k_plus_q_idx, G_idx, G′_idx)
        if G_plus_G′_idx == -1
            return 0.0
        end
        
        sum(map(1:2) do σ_idx
            c_nk_plus_q_Gσ[G_plus_G′_idx, σ_idx, :]' * c_n′k_Gσ[G′_idx, σ_idx, :]
        end)
    end)
end

"""
Some acceleration added. 
The most time consuming step here 
Since each k point has ~36000 G vectors, 
storing the indices of G+G′ in the G grid of k+q for all G′ is simply infeasible
-- also note that we haven't include the dimension for the choice of k 
and the dimension for the choice of q yet.
Since the indices of G+G' has to be re-evaluated for every (k, q, G), 
but not for every (n, n′),
we decide that calculating G+G′ for once  
and then calculating M_nn′(k, q, G) for all possible n and n′
without recalculating the positions of G+G' seems to be a good idea.
Therefore 
"""
function transition_matrix_irreducible_1BZ(wfn::BerkeleyGWSpinorWaveFunction, 
    n, n′, k_idx, q_idx, G_idx)
     
    k_plus_q_idx = find_k_plus_q_irreducible_1BZ(wfn, k_idx, q_idx)
    c_n′k_Gσ       = read_wavefunction(wfn, n′, k_idx)
    c_nk_plus_q_Gσ = read_wavefunction(wfn, n, k_plus_q_idx)
    G_plus_G′_indices = indices_of_G_plus_G′(wfn, k_idx, q_idx, G_idx)
    zero_nn′ = zeros(ComplexF64, length(n), length(n′))

    sum(map(1 : wfn.ngk[k_idx]) do G′_idx
        G_plus_G′_idx = G_plus_G′_indices[G′_idx]

        if G_plus_G′_idx == -1
            return zero_nn′ 
        end
         
        sum(map(1:2) do σ_idx
            conj(c_nk_plus_q_Gσ[G_plus_G′_idx, σ_idx, :]) * transpose(c_n′k_Gσ[G′_idx, σ_idx, :]) 
        end)
    end)
end

#endregion 
######################################################################