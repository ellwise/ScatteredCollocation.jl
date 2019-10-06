struct AffineOperator{T1<:AbstractMatrix,T2<:AbstractVector}
    A::T1
    b::T2
end

function AffineOperator(A::AbstractMatrix)
    return AffineOperator(A, zeros(eltype(A), size(A,1)))
end

(O::AffineOperator)(u) = O.A*u + Vector(O.b)
#(O::AffineOperator)(u, p, t) = O(u)
#function (O::AffineOperator)(du, u, p, t)
#    du[:] = O(u)
#    return nothing
#end
#*(O::AffineOperator, u) = O(u)
#function A_mul_B!(v, O::AffineOperator, u)
#    v[:] = O(u)
#    return nothing
#end

#Base.hcat(Os::AffineOperator...) = AffineOperator(reduce(hcat, O.A for O in Os), reduce(+, O.b for O in Os))
#Base.vcat(Os::AffineOperator...) = AffineOperator(reduce(vcat, O.A for O in Os), reduce(vcat, O.b for O in Os))
# following is not type-stable!
#function Base.hvcat(rows::Tuple{Vararg{Int}}, Os::AffineOperator...)
#    nbr = length(rows)
#    rs = Vector{Vector{AffineOperator}}(undef, nbr)
#    a = 1
#    for i = 1:nbr
#        rs[i] = reduce(hcat, Os[a:a-1+rows[i]])
#        a += rows[i]
#    end
#    return reduce(vcat, rs)
#end

# NEED TO REJIG SO THAT THE ROWS CAN BE IN AN ARBITRARY ORDER...
# THAT WAY I CAN PROJECT A BC ONTO AN ALREADY-SQUARE OPERATOR
function projbc(A::AffineOperator, B::AffineOperator, bidxs)
    m = size(A.A,1) + size(B.A,1)
    n = size(A.A,2)
    iidxs = setdiff(1:n, bidxs)
    L = spzeros(m,n)
    y = spzeros(m)
    L[iidxs,:] = A.A
    L[bidxs,:] = B.A
    y[iidxs] = A.b
    y[bidxs] = B.b
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
    C = AffineOperator(L[iidxs,iidxs], y[iidxs])
    D = AffineOperator(-L[bidxs,iidxs], -y[bidxs]) # check signs are right...
    return C, D
end

function projbc(B::AffineOperator, bidxs)
    n = size(B.A,2)
    A = AffineOperator(spzeros(n-size(B.A,1),n))
    _, D = projbc(A, B, bidxs)
    return D
end