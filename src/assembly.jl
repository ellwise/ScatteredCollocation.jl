"""
    Assembler(k, acc_order, tree) -> sca

Create a structure for assembling ScatteredCollocation differential operators.
"""
struct Assembler{T1<:Integer, T2<:NNTree}
    k::T1
    acc_order::T1
    tree::T2
end

# compute basis
# compute rhs
# compute nn
# initialise ntuple of dense output arrays (m,n where m<n)
# for each row index
#     compute an array of coeffs (for all the rhs at once)
#     for each output array
#         assign the appropriate coeffs
# rows = repeat(idxs, 1, m)'
# cols = hcat(nn...)
# map(ntuple) do vals
#     sparse(rows[:], cols[:], vals[:], n, num_points)
# return the sparse tuple
#
# eachindex(A) rather than 1:length(A)
# replace many uses of size with axes
# replace 1:length(A) with eachindex(A), or in some cases LinearIndices(A)
# replace explicit allocations like Array{Int}(size(B)) with similar(Array{Int}, axes(B))

#function derivative(T, sca::Assembler, diff_order, idxs)
#    tot_order = sum(diff_order)
#    maxorder = sca.acc_order+tot_order-1
#    # Theorem 2, Barnett's thesis
#    #unisolvent = maxorder >= (sca.k-1)/2 # assumes nodes themselves are unisolvent...
#    #if !unisolvent
#    #    error("System is not unisolvent for given derivative order. Either:
#    #           1. Increase sca.acc_order to at least $(Int(ceil((sca.k-1)/2 - tot_order)));
#    #           2. Decrease sca.k to at most $(2*(sca.acc_order+tot_order) + 1).")
#    #end
#    numpoly = countmonomials(maxorder; numdims=getdims(sca.tree))
#    stencil_size = 2*numpoly
#    xx = getpoints(sca.tree)
#    sortres = true
#    op(centre,stencil) = wlspoly(centre, stencil, diff_order; k=sca.k, npoly=numpoly)
#    nn, _ = knn(sca.tree, xx[idxs], stencil_size+numpoly, sortres)
#    nn = map(n->vcat(n[1:stencil_size-numpoly], n[end-numpoly+1:end]), nn)
#    A = assemble(T, xx, idxs, nn, op)
#    return A
#end

"""
    projbc(L, B, f, g)

This function takes an affine operator A defined on the interior
and an affine operator B defined on the boundary. It converts A 
into an equivalent operator C with no boundary dependence, whose 
solutions satisfy B. It converts B into an operator D which 
extrapolates solutions to C from the interior to the boundary.
"""
function projbc(L, y, bidxs)
    npts = size(L,1)
	iidxs = setdiff(1:npts, bidxs)
    # find a permutation to minimise subsequent fill-in
    #p = symrcm(L) # bandwidth-minimising
    p = amd(L)
    #p = colperm(L)
    #p = 1:npts
    ip = invperm(p)
    # perform Gauss-Jordan elimination, check if successful
    permute!(L, p, p)
    permute!(y, p)
    gaussjordan2!(L, y, ip[bidxs])
    permute!(L, ip, ip)
    invpermute!(y, p)
    # compute new system blocks
    # - note: A[iidxs,bidxs] is zero, A[bidxs,iidxs] is I
    A = L[iidxs,iidxs]; B = L[bidxs,iidxs]
    f = y[iidxs];       g = y[bidxs]
    return A, B, f, g
end

"""
Compute inner product weights for performing differentiation at a given point.
For multidimentional problems, each dimension can be differentiated a different
number of times.
"""
function wls(stencil, basis, rhs; weights=scaletonearest(stencil))
    
    # shift system to centre
    x = stencil .- Ref(stencil[1])
    x0 = zero(stencil[1])

    # compute the collocation problem
    W = Diagonal(weights)
    P = [g(y) for y in x, g in basis]
    p = [∂g(x0) for ∂g in rhs]

    # minimise error in weigted norm using a regularised SVD
    rtol = sqrt(eps(real(float(one(eltype(P))))))
    X = pinv(transpose(W*P), rtol=rtol)
    w = transpose(W) * X * p
	
	# ensure the weights exactly differentiate a constant
	# - weights should sum to zero, but roundoff limits this
	# - it is apparently common to adjust the central (diagonal) weight
	# - see: https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb
	# - could add a check to see how far off it was...
	# - could also look into smoothly reducing the outer weights...
	# THIS WILL CAUSE ISSUES WITH THE IDENTITY OPERATOR!
	# Maybe the identity operator shouldn't be allowed: it's just a matrix unless I add in interpolation
	# w[1] = -sum(w[2:end])
	
    return w
    
end

function scaletonearest(stencil)
    # assumes stencil is sorted
    x = stencil .- Ref(stencil[1])
    rs = norm.(x, 2)
    ε = rs[2]#rs[end]#
    w = exp.(-(rs ./ ε).^2)
end

function build(T, tree, idxs, diff_orders, accorder)

    # compute parameters
    totorder = maximum(sum.(diff_orders))
    maxorder = accorder+totorder-1
    numdims = getdims(tree)
    numpoly = countmonomials(maxorder; numdims=numdims)
    stencilsize = 2*numpoly

    # compute monomial basis and rhs
    m = MonomialExponents{numdims,Int}(numpoly)
    poly_exponents = collect(Iterators.take(m, numpoly))
    gs = Poly.(1.0, poly_exponents)
    ∂gs = [differentiate(g, ls) for g in gs, ls in diff_orders]

    # compute neighbours
    xx = getpoints(tree)
    sortres = true # find centre [1] and nearest [2] easily
    nn, _ = knn(tree, xx[idxs], stencilsize, sortres)

    # initialise ntuple of dense output arrays (m,n where m<n)
    denseops = ntuple(length(diff_orders)) do x
        Array{T,2}(undef, stencilsize, length(idxs))
    end

    # compute dense arrays of op values
    progress = 0
    print("Assembling...")
    for j in eachindex(idxs)
        stencil = xx[nn[j]]
        coeffs = wls(stencil, gs, ∂gs)
        for i in eachindex(denseops)
            denseops[i][:,j] = coeffs[:,i]
        end
        progress += 1
        if progress > length(idxs)/10.0
            print("+")
            progress = 0
        end
    end
    println("*")

    # convert to sparse arrays for output
    rows = repeat((1:length(idxs))', stencilsize, 1)
    cols = reduce(hcat, nn) # hcat(nn...) -> splat is not good for compilation
    sparseops = map(denseops) do vals
        sparse(rows[:], cols[:], vals[:], length(idxs), length(xx))
    end
    return sparseops
    #
    # replace many uses of size with axes
    # replace 1:length(A) with eachindex(A), or in some cases LinearIndices(A)
    # replace explicit allocations like Array{Int}(size(B)) with similar(Array{Int}, axes(B))
end
# Wishlist:
# - CLArray{T}
# - CuArray{T}
# - Array{T}?

"""
struct AffineOp{T1,T2} where {T1<:AbstractMatrix, T2<:AbstractVector}
    A::T1
    b::T2
end

(L::AffineOp)(u, p, t) = L.A*u + L.b
function (L::AffineOp)(du, u, p, t)
    du[:] = L(u, p, t)
end
# update_coefficients!(L,u,p,t)
is_constant(L::AffineOp) = true
"""