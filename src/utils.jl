
"""
    Compute the 2-norm.
"""
norm2(x) = sqrt(sum(x.^2))

"""
    Pochhammer symbol (rising factorial).
"""
pochhammer(x::Integer, n::Integer) = n==0 ? one(x) : (x+n-1)*pochhammer(x, n-1)

# copy tuples with modifications
inctuple(tup, idx, val) = typeof(tup)(idx==j ? tup[j] + val : tup[j] for j=1:length(tup))
settuple(tup, idx, val) = inctuple(tup, idx, val - tup[idx])
# typeof(tup)(idx==j ? val : tup[j] for j=1:length(tup))

# make this an iterator!
function halton_sequence(i::Integer, b::Integer)
    f = 1.0
    r = 0.0
    while i > 0
        f /= b
        r += f*mod(i, b)
        i = convert(Int, fld(i, b))
    end
    return r
end

# https://people.sc.fsu.edu/~jburkardt/m_src/monomial/monomial.html
"""
    unique_monomials(ndims, maxorder)

Return the number of unique monomials up to a given order.

# Examples
```jldoctest
julia> unique_monomials(2, 3)
10
```
"""
unique_monomials(numdims, maxorder) = binomial(maxorder + numdims, maxorder)
#function unique_monomials(numdims::Int, maxorder::Int)
#    @assert numdims > 0
#    @assert maxorder ≥ 0
#    binomial(maxorder + numdims, maxorder)
#end

"""
    monomial_exponents(ndims)

Return an iterable that generates monomial exponents.
"""
function monomial_exponents(numdims::Int)
    @assert numdims > 0
    if numdims==1
        return IterTools.iterated(grlex, 0)
    end
    return IterTools.iterated(grlex, ntuple(x->zero(x), numdims))#zeros(Int, numdims))
end

"""
Return the next set of monomial exponents given a current set. Uses graded
lexicographic ordering.
"""
grlex(x_in::Integer) = x_in + 1
function grlex(x_in::Tuple{Vararg{Integer}})#Vector{Int})
    # be careful with modifying arrays! In this case, I wanted to return a copy...
    #x = copy(x_in) # deepcopy not needed because eltype is immutable
    x = x_in
    m = length(x)
    @assert m ≥ 0
    #@assert all(x .≥ 0) # ensure that 0 <= x[i]
    # find the index of the rightmost nonzero entry of X.
    i = findlast(y->y≠0, x)
    # set T = X(I), set X(I) to zero, increase X(I-1) by 1, increment X(M) by T-1.
    if i == nothing
        x = settuple(x, m, 1)#x[m] = 1
        return x
    elseif i == 1
        t = x[1] + 1
        im1 = m
    elseif i > 1
        t = x[i]
        im1 = i - 1
    end
    x = settuple(x, i, 0)#x[i] = 0
    x = inctuple(x, im1, 1)#x[im1] +=  1
    x = inctuple(x, m, t-1)#x[m] += t - 1
    return x
end

"""
Return the next set of monomial exponents given a current set. Uses graded
_reverse_ lexicographic ordering.
"""
grevlex(x_in::Integer) = x_in + 1
function grevlex(x_in::Tuple{Vararg{Integer}})
    # be careful with modifying arrays! In this case, I wanted to return a copy...
    x = x_in#copy(x_in) # deepcopy not needed because eltype is immutable
    m = length(x)
    @assert m ≥ 0
    # ensure that exponents are equal to or beyond the initial grevlex term
    @assert all(x .≥ 0)
    # Seek the first index 1 < I for which 0 < X(I).
    j = findfirst(y->y>0, x[2:m])#view(x, 2:m))
    j = j==nothing ? 1 : j + 1 # correct the indexing
    if j == 1
        t = x[1]
        x = settuple(x, 1, 0)#x[1] = 0
        x = settuple(x, m, t+1)#x[m] = t + 1
    elseif j < m
        x = inctuple(x, j, -1)#x[j] -= 1
        t = x[1] + 1
        x = settuple(x, 1, 0)#x[1] = 0
        x = inctuple(x, j-1, t)#x[j-1] += t
    elseif j == m
        t = x[1]
        x = settuple(x, 1, 0)#x[1] = 0
        x = settuple(x, j-1, t+1)#x[j-1] = t + 1
        x = inctuple(x, j, -1)#x[j] -= 1
    end
    return x
end
