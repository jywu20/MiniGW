######################################################################
#region IO

export Ry_BGW, 
    BerkeleyGWSpinorWaveFunction,
    read_wavefunction,
    read_wavefunctions, 
    transition_matrix_def

const Ry_BGW =  13.6056925

struct BerkeleyGWSpinorWaveFunction <: AbstractGWWaveFunction
    fid::HDF5.File
    
    nrk::Int
    irreducible_1BZ::Matrix{Float64}
    full_1BZ::Matrix{Float64}
    el::Matrix{Float64}

    gk_ranges::Vector{UnitRange} 
    ngk::Vector{Int} 
    ngkmax::Int
    gvecs::Vector{Matrix{Int}}
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
    full_1BZ = expand_sym_reduced_grid(rk, map(1 : ntran) do sym_idx
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
    
    gvecs_original_form = fid["wfns/gvecs"] |> read
    gvecs = Vector{Matrix{Float64}}(undef, nrk)
    for k_idx in 1 : nrk
        gvecs[k_idx] = gvecs_original_form[:, gk_ranges[k_idx]]
    end
    
    BerkeleyGWSpinorWaveFunction(fid, 
        nrk, rk, full_1BZ, 
        el, 
        gk_ranges, ngk, ngkmax, gvecs)
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
#region Transition matrix 

"""
`q_idx` is assumed to be the index in `wfn.irreducible`; 
some optimization may be done here.

TODO: this is actually wrong - 
at this or that point, 
we will need to consider the possibility that 
k and q are in the irreducible 1BZ 
but k + q is not.
"""
function find_k_plus_q(wfn::BerkeleyGWSpinorWaveFunction, k_idx, q_idx)
    irreducible_1BZ = wfn.irreducible_1BZ
    k_plus_q =  irreducible_1BZ[:, k_idx] + irreducible_1BZ[:, q_idx]
    find_in_periodic_grid(irreducible_1BZ, k_plus_q)
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
function find_G_plus_G′(wfn::BerkeleyGWSpinorWaveFunction, k_idx, k_plus_q_idx, G_idx, G′_idx)
    G_grid_of_k = wfn.gvecs[k_idx] 
    G_grid_of_k_plus_q = wfn.gvecs[k_plus_q_idx]

    G  = G_grid_of_k[:, G_idx]
    G′ = G_grid_of_k[:, G′_idx]
    find_in_grid(G_grid_of_k_plus_q, G + G′)
end

"""
The most naive implementation of the transition matrix M_nn'(k, q, G).

TODO: testing 
"""
function transition_matrix_def(wfn::BerkeleyGWSpinorWaveFunction, 
    n, n′, k_idx, q_idx, G_idx)
    
    k_plus_q_idx = find_k_plus_q(wfn, k_idx, q_idx)
    
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

#endregion 
######################################################################