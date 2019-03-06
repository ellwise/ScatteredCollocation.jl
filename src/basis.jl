# fractional derivatives: https://arxiv.org/pdf/1612.07563.pdf
# https://en.wikipedia.org/wiki/Fractional_calculus#Fractional_derivative_of_a_basic_power_function

"""
    poly(x, exponents)
    poly(x, exponents, diff_orders)

Compute a specified monomial basis function or its derivatives. Multidimensional
arguments should be given as tuples.
"""
function poly end
poly(x, exponent) = x^exponent
function poly(x, exponent, diff_order)
    if diff_order ≤ exponent
        # https://www.wolframalpha.com/input/?i=d%5Ek%2Fdx%5Ek+x%5Ea
        new_exponent = exponent-diff_order
        coeff = pochhammer(1+new_exponent, diff_order)
        return coeff * poly(x, new_exponent)
    end
    return zero(x)
end
poly(x::Tuple, exponents::Tuple) = prod(poly.(x, exponents)) # @assert all(exponents .≥ 0)
poly(x::Tuple, exponents::Tuple, diff_orders::Tuple) = prod(poly.(x, exponents, diff_orders))

"""
    phs(x, k)
    phs(x, k, exponents)
    phs(x, k, exponents, diff_orders)

Compute an odd polyharmonic spline (PHS) or its derivatives. If exponents are
given, the PHS is defined by poly(x, exponents) * phs(x, k). This facilitates
a recursive definition of the PHS's derivatives. Multidimensional arguments
should be given as tuples.
"""
function phs end
function phs(x, k)
    # to ensure the system is uniquely solvable, Barnet gives theorems for only r^k phs
    # I don't need to worry about replicating terms in 1D, because I always choose a maximal k, i.e. k > maxorder
    #@assert isodd(k)
    #@assert k > 0
    return abs(x)^k
end
phs(x, k, exponent) = poly(x, exponent) * phs(x, k)
function phs(x, k, exponent, diff_order)
    if diff_order > 0
        if exponent ≤ 0
            left_term = zero(x)
        else
            left_term = exponent*phs(x, k, exponent-1, diff_order-1)
        end
        right_term = k*phs(x, k-2, exponent+1, diff_order-1)
        return left_term + right_term
    end
    return phs(x, k, exponent)
end

phs(x::Tuple, k) = phs(norm2(x), k)
phs(x::Tuple, k, exponents::Tuple) = poly(x, exponents) * phs(x, k)
function phs(x::Tuple, k, exponents::Tuple, diff_orders::Tuple)
    # add a diff_orders check to ensure the PHS is high enough order?
    # OR, I NEED TO REDEFINE phs SO THAT THIS ISN'T A PROBLEM ANY MORE
    dim = findfirst(y->y≠0, diff_orders)
    if dim≠nothing
        diff_orders_new = inctuple(diff_orders, dim, -1)#copy(diff_orders)
        #diff_orders_new[dim] -= 1
        coeff = exponents[dim]
        if coeff ≤ 0
            left_term = zero(eltype(x))
        else
            exponents_new = inctuple(exponents, dim, -1)#copy(exponents)
            #exponents_new[dim] -= 1
            left_term = coeff * phs(x, k, exponents_new, diff_orders_new)
        end
        # do I need to deal with k=0 specially?
        exponents_new = inctuple(exponents, dim, 1)#copy(exponents)
        #exponents_new[dim] += 1
        right_term = k * phs(x, k-2, exponents_new, diff_orders_new)
        return left_term + right_term
    end
    return phs(x, k, exponents)
end
