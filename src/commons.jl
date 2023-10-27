export find_in_grid, wigner_seitz, find_in_periodic_grid, 
    equivalent_grids,
    equivalent_vectors,
    expand_sym_reduced_grid

"""
Find `vec` in `grid`. The return value is the position of `vec` in `grid`.
`grid[:, 1]` is the first vector it contains, and so on.
The index of the second dimension of `grid` is assumed to start from 1.
When `vec` can't be found, -1 is returned.
"""
function find_in_grid(grid::AbstractMatrix, vec::AbstractVector; tol = 1e-10)
    n_vecs = size(grid)[2]
    
    if size(grid)[1] != length(vec)
        throw("Dimension mismatch: the grid contains $(size(grid)[1])-dimensional vectors while vec is $(length(vec))-dimensional.")
    end
    
    for idx_vec in 1 : n_vecs
        if norm(grid[:, idx_vec] - vec) <= tol  
            return idx_vec
        end
    end

    return -1
end

function find_in_grid(grid::Matrix{I}, vec::Vector{I}) where {I <: Integer}
    n_vecs = size(grid)[2]
    
    if size(grid)[1] != length(vec)
        throw("Dimension mismatch: the grid contains $(size(grid)[1])-dimensional vectors while vec is $(length(vec))-dimensional.")
    end
    
    for idx_vec in 1 : n_vecs
        if grid[:, idx_vec] == vec
            return idx_vec
        end
    end

    return -1
end

"""
Move a number to the "1BZ" near zero, 
assuming that x is equivalent to x + 1.
-0.5 is moved to 0.5.
"""
wigner_seitz(x::Number) = x - ceil(x - 0.5)

"""
Find `vec` in `grid`, but assuming that x and x + 1 are equivalent. 
The return value is the position of `vec` in `grid`.
`grid[:, 1]` is the first vector it contains, and so on.
The index of the second dimension of `grid` is assumed to start from 1.
When `vec` can't be found, -1 is returned.
"""
function find_in_periodic_grid(grid::AbstractMatrix, vec::AbstractVector; tol = 1e-10)
    n_vecs = size(grid)[2]
    
    if size(grid)[1] != length(vec)
        throw("Dimension mismatch: the grid contains $(size(grid)[1])-dimensional vectors while vec is $(length(vec))-dimensional.")
    end
    
    for idx_vec in 1 : n_vecs
        if norm(wigner_seitz.(grid[:, idx_vec] - vec)) <= tol  
            return idx_vec
        end
    end

    return -1
end

function equivalent_grids(grid1::AbstractMatrix, grid2::AbstractMatrix; tol = 1e-10)
    if size(grid1) != size(grid2)
        return false
    end
    
    n_vecs = size(grid1)[2]

    for vec_idx in 1 : n_vecs
        if find_in_grid(grid1, grid2[:, vec_idx], tol = tol) == -1
            return false
        end
        
        if find_in_grid(grid2, grid1[:, vec_idx], tol = tol) == -1
            return false
        end
        
        return true
    end
end

"""
Expand `grid` according to the list of symmetry operations `operations`.
`grid[:, 1]` is the first vector it contains, and so on.
The index of the second dimension of `grid` is assumed to start from 1.

If a point obtained from a symmetry operation is decided to be the same as another point 
then the former is not included.

This function is designed for 1BZ expansion so things like fractional translation 
    are not considered (yet).
"""
function expand_sym_reduced_grid(grid::AbstractMatrix, operations::AbstractVector; tol = 1e-10)
    n_vecs = size(grid)[2]
    results = AbstractVector[]
    full_to_irreducible = Int[] 
    
    for idx_k_irreducible in 1 : n_vecs
        current_vec = grid[:, idx_k_irreducible]
        push!(results, current_vec)
        push!(full_to_irreducible, idx_k_irreducible)

        for op in operations
            vec_after_sym = op * current_vec
            
            found_duplicate = false
            for prev_vec in results
                if norm(wigner_seitz.(prev_vec - vec_after_sym)) <= tol
                    found_duplicate = true 
                    break
                end
            end
            
            if ! found_duplicate
                push!(results, vec_after_sym)
                push!(full_to_irreducible, idx_k_irreducible)
            end
        end
    end
    
    hcat(results...), full_to_irreducible
end

function equivalent_vectors(v1::AbstractVector, v2::AbstractVector; tol = 1e-10)
    norm(v1 - v2) <= tol
end

function equivalent_periodic_vectors(v1::AbstractVector, v2::AbstractVector; tol = 1e-10)
    norm(wigner_seitz.(v1 - v2)) <= tol  
end

#region Wave function analysis tools 

export top_n_positions

function top_n_positions(vec::AbstractVector, n::Integer)
    partialsortperm(abs.(vec), 1 : n, rev = true)
end

function main_components(ψ::AbstractArray, n::Integer)
    ψ[top_n_positions(ψ, n)]
end

#endregion 