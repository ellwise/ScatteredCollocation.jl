
# Wishlist:
# - Array{T}
# - SparseMatrixCSC{Ti, Tv}
# - CLArray{T}
# - CuArray{T}

function assemble(::Type{Array{T}}, x, ind, nn, op) where T<:AbstractFloat
    xx = typeof(x)==T ? x : map(y->convert.(T, y), x)
    num_points = length(x)
    num_dof = length(nn)
    L = zeros(Float64, num_dof, num_points)
    for row = 1:num_dof
        L[row, nn[row]] = op(xx[ind[row]], xx[nn[row]])
    end
    return L
end

function assemble(::Type{SparseMatrixCSC{T}}, x, ind, nn, op) where T<:AbstractFloat
    xx = typeof(x)==T ? x : map(y->convert.(T, y), x)
    num_points = length(x)
    num_dof = length(nn)
    num_nz = sum(length.(nn))
    ii = Vector{Int}(undef, num_nz)
    jj = Vector{Int}(undef, num_nz)
    vv = Vector{Float64}(undef, num_nz)
    k = 0
    for row = 1:num_dof
        stencil_size = length(nn[row])
        ii[k+1:k+stencil_size] .= row
        jj[k+1:k+stencil_size] = nn[row]
        vv[k+1:k+stencil_size] = op(xx[ind[row]], xx[nn[row]])
        k += stencil_size
    end
    L = sparse(ii, jj, vv, num_dof, num_points)
    return L
end
