"""
    Poly(coeff, exponents)

Defines a monomial with an associated coefficient. Once constructed, the Poly
can be evaluated like a function, or differentiated to yield another Poly.
"""
struct Poly
    coeff::Float64
    exponents::Vector{Int}
end

# not strictly right (doesn't deal with 1/x correctly), but helps with differentiation
function (f::Poly)(x::AbstractVector)
    val = f.coeff * prod(x.^f.exponents)
    return iszero(f.coeff) ? zero(val) : val
    # below is not type-stable (eltype(x)) not neccessarily typeof(f.coeff * prod(x.^f.exponents))
    #return iszero(f.coeff) ? zero(eltype(x)) : f.coeff * prod(x.^f.exponents)
end

function differentiate(f::Poly, diff_orders)
    dim = findfirst(y->y≠0, diff_orders)
    if isnothing(dim)
        return f
    else
        coeff = f.coeff
        ps = f.exponents
        diff_orders_new = copy(diff_orders); diff_orders_new[dim] -= 1
        psm1 = copy(ps); psm1[dim] -= 1
        g = Poly(coeff*ps[dim], psm1)
        return differentiate(g, diff_orders_new)
    end
end

"""
    PHS(k, islog)

Defines a generalised polyharmonic spline. Once constructed, the PHS
can be evaluated like a function, or differentiated to yield a Vector
of PHSPoly.
"""
struct PHS
    k::Int
    islog::Bool
end

function (f::PHS)(x::AbstractVector)
    r = norm(x, 2)
    if f.islog
        return r < 1 ? r^(f.k-1)*log(r^r) : r^f.k*log(r)
    else
        return r^f.k
    end
end

function differentiate(f::PHS, diff_orders)
    g = Poly(1.0, zero(diff_orders))
    h = PHSPoly(f, g)
    return differentiate(h, diff_orders)
end

"""
    PHSPoly(f, g)

Defines a product of a monomial and a polyharmonic spline. Once constructed, the 
PHSPoly can be evaluated like a function, or differentiated to yield a Vector of
PHSPoly.

This basis function is useful because it is closed under differentiation. This
feature is used to form a recursive differentiation algorithm.
"""
struct PHSPoly
    phs::PHS
    poly::Poly
end

function (f::PHSPoly)(x::AbstractVector)
    val = f.phs(x) * f.poly(x)
    # most NaN/Inf base cases have been avoided during the recursion
    return isnan(val) ? zero(val) : val
end

(fs::Vector{PHSPoly})(x::AbstractVector) = mapreduce(f->f(x), +, fs)

function differentiate(f::PHSPoly, diff_orders)
    dim = findfirst(y->y≠0, diff_orders)
    if isnothing(dim)
        return [f]
    else
        coeff = f.poly.coeff
        ps = f.poly.exponents
        k = f.phs.k
        diff_orders_new = copy(diff_orders); diff_orders_new[dim] -= 1
        psm1 = copy(ps); psm1[dim] -= 1
        psp1 = copy(ps); psp1[dim] += 1
        # the iszero(...) checks below prevent all of the NaN base cases
        if f.phs.islog
            term1 = iszero(ps[dim]) ? PHSPoly[] : differentiate(PHSPoly(PHS(k, true),    Poly(coeff*ps[dim], psm1)), diff_orders_new) 
            term2 = iszero(k)       ? PHSPoly[] : differentiate(PHSPoly(PHS(k-2, true),  Poly(coeff*k,       psp1)), diff_orders_new)
            term3 =                               differentiate(PHSPoly(PHS(k-2, false), Poly(coeff,         psp1)), diff_orders_new)
            return vcat(term1, term2, term3)
        else
            term1 = iszero(ps[dim]) ? PHSPoly[] : differentiate(PHSPoly(PHS(k, false),   Poly(coeff*ps[dim], psm1)), diff_orders_new)
            term2 = iszero(k)       ? PHSPoly[] : differentiate(PHSPoly(PHS(k-2, false), Poly(coeff*k,       psp1)), diff_orders_new)
            return vcat(term1, term2)
        end
    end
end