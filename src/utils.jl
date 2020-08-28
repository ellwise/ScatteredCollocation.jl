"""
    p = colperm(A)

Find a permutation that sorts columns in a sparse matrix by their non-zero count.
"""
function colperm(A::SparseMatrixCSC)
    colcounts = [A.colptr[j+1]-A.colptr[j] for j in 1:A.n]
    p = sortperm(colcounts)
    return p
end

"""
    n = fnz(A)

Compute the fraction of non-zero elements in a sparse array.
"""
fnz(A::SparseMatrixCSC) = nnz(A)/prod(size(A))

"""
    scplot!(scene, tree, u; kwargs...)

Plot a model solution defined on scattered points using Makie.
"""
function scplot!(scene, tree, u; kwargs...)
    numdims = getdims(tree)
    x = getpoints(tree)
    sortres = true
    _, ds = knn(tree, getpoints(tree), 2, sortres)
    Δxs = [d[2] for d in ds]
    if numdims==1
        lines!(scene, map(y->y[1], x), u; kwargs...)
    elseif numdims==2
        extrakwargs = Dict(:markersize=>1.1*Δxs, :color=>u)
        scatter!(scene, x; extrakwargs..., kwargs...)
    #elseif numdims==3
    #    extrakwargs = Dict(:markersize=>1.1*Δxs, :color=>u, :alpha=>0.99)#abs.(u))
    #    scatter!(scene, x; extrakwargs..., kwargs...)
    end
end

"""
    scspy!(scene, A)

Generate a spy plot for a sparse matrix using Makie.
"""
function scspy!(scene, A::SparseMatrixCSC)
    B = zeros(Int, size(A))
    rows = rowvals(A)
    ncols = size(A,2)
    for col in 1:ncols
        for j in nzrange(A, col)
            row = rows[j]
            B[row,col] = 1
        end
    end
    heatmap!(scene, transpose(B)[:,end:-1:1])
    #clims = maximum(abs.(A))*[-1, 1]
    #heatmap!(scene, transpose(Array(A))[:,end:-1:1]; colorrange=clims, colormap=:pu_or)
end

"""
    haltonsequence(i, b; f=1.0, r=0.0)

Compute low-discrepancy Halton sequences.
"""
function haltonsequence(i::Integer, b::Integer; f=1.0, r=0.0)
    if iszero(i)
        return r
    else
        f /= b
        r += f*mod(i, b)
        i = convert(Int, fld(i, b))
        return haltonsequence(i, b, f=f, r=r)
    end
end

# https://people.sc.fsu.edu/~jburkardt/m_src/monomial/monomial.html
"""
    countmonomials(maxorder; ndims)

Return the number of unique monomials up to a given order.

# Examples
```jldoctest
julia> unique_monomials(2, 3)
10
```
"""
countmonomials(maxorder; numdims) = binomial(maxorder + numdims, maxorder)

"""
    MonomialExponents{N, T}(npoly)

This is an iterable that generates monomial exponents. N and T are Integers
that specify the number of dimensions, and the type of exponent values. The
argument _npoly_ specifies the maximum number of terms that the iterable returns.
"""
struct MonomialExponents{N,T}
    npoly::Integer
end
Base.iterate(iter::MonomialExponents{N,T}) where {N,T} = ntuple(x->zeros(MVector{N,T}), 2)
Base.iterate(iter::MonomialExponents{N,T}, exponents::MVector{N,T}) where {N,T} = ntuple(x->next_exponent(exponents), 2)
Base.length(m::MonomialExponents) = m.npoly

"""
Return the next set of monomial exponents given a current set. Uses either
graded or _reverse_ graded lexicographic ordering (ordering = :grlex/:grevlex).
"""
function next_exponent(x_in::AbstractVector; ordering=:grlex)
    # be careful with modifying arrays! In this case, I wanted to return a copy...
    x = copy(x_in) # deepcopy not needed because eltype is immutable
    m = length(x)
    @assert m ≥ 0
    # ensure that exponents are equal to or beyond the initial grevlex term
    @assert all(x .≥ 0)
    if ordering==:grlex
        # find the index of the rightmost nonzero entry of X.
        i = findlast(y->y≠0, x)
        # set T = X(I), set X(I) to zero, increase X(I-1) by 1, increment X(M) by T-1.
        if isnothing(i)
            x[m] = 1
            return x
        elseif i == 1
            t = x[1] + 1
            im1 = m
        elseif i > 1
            t = x[i]
            im1 = i - 1
        end
        x[i] = 0
        x[im1] +=  1
        x[m] += t - 1
    elseif ordering==:grevlex
        # Seek the first index 1 < I for which 0 < X(I).
        j = findfirst(y->y>0, view(x, 2:m))
        j = isnothing(j) ? 1 : j + 1 # correct the indexing
        if j == 1
            t = x[1]
            x[1] = 0
            x[m] = t + 1
        elseif j < m
            x[j] -= 1
            t = x[1] + 1
            x[1] = 0
            x[j-1] += t
        elseif j == m
            t = x[1]
            x[1] = 0
            x[j-1] = t + 1
            x[j] -= 1
        end
    end
    return x
end