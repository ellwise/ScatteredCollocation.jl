"""
    gaussjordan(A, b, idxs)

Compute a Gauss-Jordan elimination matrix for a given set of column indices.
"""
function gaussjordan!(A::SparseMatrixCSC, y::AbstractVector, idxs)
    progress = 0
    nidxs = length(idxs)
    print("Performing Gauss-Jordan elimination...")
    # be careful! I need idxs to be in order, so I'm always looking below the main diagonal
    while ~isempty(idxs)
		# find the pivot indices (k,j)
        #pivoting = :partial
		#if pivoting==:nopivot
			j = idxs[1]
			k = j
		#elseif pivoting==:partial
		#	j = idxs[1]
		#	scales = maximum(abs.(A[idxs,:]); dims=2)
		#	k = idxs[argmax(abs.(A[idxs,j]) ./ scales)]
		#end
        # divide the pivot row by the pivot element
        p = A[k,j]
		A[k,:] ./= p
		y[k:k] ./= p
        # subtract multiples of the pivot row from all other rows
        # note: an element-by-element algorithm is slow because sparse matrices are immutable (hence, each change allocates a new matrix)
        prow = A[k:k,:]
        pcol = A[:,j]
        U = pcol .* prow
        u = y .* y[k]
        U[k,:] .= 0 # don't update the pivot row
        u[k] = 0
        #rtol = sqrt(eps(real(float(one(eltype(A))))))
        reducefillin!(U)#; byrow=true)#, rtol=rtol)
        A[:,:] -= U
		y[:] -= u
        # indicate progress
        progress += 1
        if progress > nidxs/10.0
            print("+")
            progress = 0
        end
        idxs = idxs[2:end]#filter!(x->x≠j, idxs) FILTER VERSION IS NOT TYPE STABLE
    end
    println("*")
    # check residuals
    nidxs = setdiff(1:size(A,1), idxs)
    R1 = A[nidxs,idxs]
    R2 = A[idxs,idxs]
    if ~(R1≈zero(R1) && R2≈I)
        @warn "Failed to transform specified columns to reduced echelon form."
    end
    return nothing
end

function gaussjordan2!(L, y, bidxs; tol=1e-6)
	# extract input blocks
    npts = size(L,1)
    nbnd = length(bidxs)
    iidxs = setdiff(1:npts, bidxs)
	Ab = L[iidxs,bidxs]
    Bb = L[bidxs,bidxs]
    
    LHS = sparse(Bb')
    RHS = spzeros(nbnd, npts)
    RHS[:,iidxs] = -Ab'
    RHS[:,bidxs] = sparse(1.0I, nbnd, nbnd)
    LU = ilu(LHS, τ=tol)
    # compute Gauss-Jordan columns
    progress = 0.0
    print("Performing Gauss-Jordan elimination...")
	Mrows = map(1:npts) do j
        Mj = gmres(LHS, Vector(RHS[:,j]); Pl=LU, tol=tol)
        # indicate progress
        progress += 1.0
        if progress > npts/10.0
            print("+")
            progress = 0
        end
		return reducefillin!(sparse(Mj))'
	end
	M = reduce(vcat, Mrows)
    println("*")
	# assemble and apply Gauss-Jordan matrix
    GJ = sparse(1.0I, npts, npts)
	GJ[:,bidxs] = M
    L[:,:] = GJ*L
    reducefillin!(L)
	y[:] = GJ*y
	# check residuals
    R1 = L[iidxs,bidxs]
    R2 = L[bidxs,bidxs]
    if ~(isapprox(R1, zero(R1); atol=tol) && isapprox(R2, I; atol=tol, rtol=tol))
        @warn "Failed to transform specified columns to reduced echelon form."
    end
	return nothing
end


"""
    reducefillin!(A; byrow=false, atol=0.0, rtol::Real=atol>0 ? 0 : n*ϵ)

Drop small values from a SparseMatrixCSC. See pinv for advice around tolerances.
"""
function reducefillin!(A::SparseMatrixCSC{T}; byrow::Bool=false, atol::Real=0.0, rtol::Real=(eps(real(float(one(T))))*min(size(A)...))*iszero(atol)) where T
    tol =  byrow ? max.(rtol*maximum(A;dims=2), atol) : max(rtol*maximum(A), atol)
    rows = rowvals(A)
    vals = nonzeros(A)
    for col = 1:size(A,2)
        for j in nzrange(A,col)
            row = rows[j]
            val = vals[j]
            smallval = byrow ? abs(val)<tol[row] : abs(val)<tol
            A[row,col] = smallval ? zero(val) : val
        end
    end
    dropzeros!(A)
end

function reducefillin!(v::SparseVector{T}; atol::Real=0.0, rtol::Real=(eps(real(float(one(T))))*min(size(v)...))*iszero(atol)) where T
    tol =  max(rtol*maximum(v), atol)
    for j in v.nzind
        if abs(v[j]) < tol
            v[j] = 0
        end
    end
    dropzeros!(v)
end

#### from Matlab: rref (reduced row echelon form)
#### Pivot tolerance, specified as a scalar.
#### If the largest element (by absolute value) in a pivot column is below the tolerance, then the column is zeroed out.
#### This prevents division and multiplication with nonzero pivot elements smaller than the tolerance.
#### max(size(A))*eps*norm(A,inf) (default)
