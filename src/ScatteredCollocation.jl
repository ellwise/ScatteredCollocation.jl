"""
    Scattered collocation provides methods for generating inner product weights
    corresponding to differential operators on scattered stencil points.
"""
module ScatteredCollocation

import IterTools
using SparseArrays

include("./geometry.jl")
include("./basis.jl")
include("./utils.jl")
include("./assembly.jl")
include("./dolfin.jl")

export differentiate, identityop, gradientop, laplacianop, linearop, assemble

"""
Identity operator. Will perform interpolation if the centre node is not in the
set of stencil nodes.
"""
identityop(centre, stencil) = differentiate(centre, stencil, floor.(Int, zero.(centre)))

"""
Laplacian operator.
"""
laplacianop(centre::Number, stencil) = differentiate(centre, stencil, 2)
function laplacianop(centre, stencil)
    ndims = length(centre)
    diff_orders = [settuple(floor.(Int, zero.(centre)), j, 2) for j=1:ndims]
    d2tup = differentiate(centre, stencil, diff_orders)
    return sum(d2tup)
end

"""
Gradient operator. Will return inner product weights for each dimension.
"""
gradientop(centre::Number, stencil) = differentiate(centre, stencil, 1)
function gradientop(centre, stencil)
    ndims = length(centre)
    diff_orders = [settuple(floor.(Int, zero.(centre)), j, 1) for j=1:ndims]
    d1tup = differentiate(centre, stencil, diff_orders)
    return d1tup
end

"""
General linear operator. Returns a linear combination of differential operators.
"""
function linearop(centre, stencil, diff_orders, coefficients)
    weights = differentiate(centre, stencil, diff_orders)
    return sum(coefficients.*weights)
end

"""
Compute inner product weights for performing differentiation at a given point.
For multidimentional problems, each dimension can be differentiated a different
number of times.
"""
function differentiate(centre::T, stencil::Vector{T}, diff_orders) where T

    # check whether order and length(stencil) are compatible
    # @assert order < length(stencil)
    # ADD A CHECK FOR COINCIDENT NODES! MAYBE HERE, MAYBE HIGHER. MY TESTS WERE FAILING BECAUSE OF THIS!
    # ADJUST MINNODES TO ACCOUNT FOR DERIVATIVE ORDRR!

    # validate inputs
    numdims, num_nodes = _checkinputs(centre, stencil, diff_orders)

    # shift the whole system to the origin
    x0 = zero.(centre)
    x = reshape(broadcast((a,b)->a.-b, copy(stencil), Ref(centre)), :, 1)
    #x = reshape([s.-centre for s in stencil], :, 1)
    #x = reshape(map(.-, copy(stencil), Ref(centre)), :, 1)

    # determine an appropriate set of basis functions
    k, poly_exponents = _choosebasis(numdims, num_nodes)

    # compute the blocks of the collocation problem
    A, P = _lhs(x, k, poly_exponents)
    a, p = _rhs(x0, x, k, poly_exponents, diff_orders)

    # concatenate into asingle least-squares system
    num_poly = size(p, 1)
    num_nodes = length(stencil)
    Z = zeros(eltype(p), num_poly, num_poly)
    B = [A P; P' Z]
    f = [a; p]

    # solve and retain the appropriate coefficients
    c = B \ f
    c = size(c, 2) == 1 ? c[1:num_nodes] : tuple([c[1:num_nodes, j] for j = 1:size(c, 2)]...)
    return c

end

"""
Determine the smallest feasible stencil for computing differentiation weights.
"""
function smalleststencil(numdims::Int, total_difforder::Int)
    # The PHS order must be at least 3 for differentiation (as defined below) to
    # work. For unisolvency? we then need polynomial terms up to order l >= (k-1)/2.
    # Also need twice as many nodes as monomial terms
    # This leads to a minimum stencil size.
    # NOT SURE THAT THIS RULE IS TRUE FOR DIFFERENTIATION, ONLY JUSTIFIED FOR INTERPOLATION!
    # need to make sure the PHS is high-enough order for differentiation... hence total_difforder
    min_phs_order = 3 + 2*total_difforder
    min_maxorder = convert(Int, cld(min_phs_order-1, 2))
    min_num_nodes = 2 * unique_monomials(numdims, min_maxorder)
end

function _checkinputs(centre, stencil, diff_orders)

    # The PHS order must be at least 3 for differentiation (as defined below) to
    # work. For unisolvency? we then need polynomial terms up to order (k-1)/2.
    # This leads to a minimum stencil size.
    numdims = length(centre)
    num_nodes = length(stencil)
    #min_num_nodes = smalleststencil(numdims, sum(diff_orders)) # THIS DOESN'T WORK ANYMORE! WAS DESIGNED AROUND DIFF_ORDERS HAVING JUST ONE!
    #@assert num_nodes ≥ min_num_nodes
    return numdims, num_nodes

end

function _choosebasis(numdims, num_nodes)

    # Now, the PHS order should be as high as possible, given the requirements
    # of unisolvency. Specifically, the highest order for which all monomials
    # are included is found, and the PHS order is based on this
    num_poly = convert(Int, fld(num_nodes, 2))
    maxorder = sum(IterTools.nth(monomial_exponents(numdims), num_poly))
    if num_poly < unique_monomials(numdims, maxorder)
        maxorder -= 1
    end
    k = 2*maxorder + 1
    # BARNETT'S THESIS ONLY CONSIDERS ||x||^m, NOT ALL PHS! However, the formula above already ensures m is odd

    poly_exponents = reshape(collect(Iterators.take(monomial_exponents(numdims), num_poly)), :, 1)

    return k, poly_exponents

end

function _lhs(x, k, poly_exponents)

    # compute the PHS interpolation matrix and monomial matrix
    # SHOULD LEVERAGE SYMMETRY HERE! Distances.jl uses for loops and in-place operations
    A = phs.(broadcast((a,b)->a.-b, x, reshape(x, 1, :)), k)
    P = poly.(x, reshape(poly_exponents, 1, :))
    return A, P

end

function _rhs(x0, x, k, poly_exponents, diff_orders::Vector)
    phs_exponents = zero.(x0) # deals with arrays and scalars
    a = phs.(broadcast((a,b)->a.-b, Ref(x0), x), k, Ref(phs_exponents), reshape(diff_orders, 1, :))
    p = poly.(Ref(x0), poly_exponents, reshape(diff_orders, 1, :))
    return a, p
end

function _rhs(x0, x, k, poly_exponents, diff_orders)
    phs_exponents = zero.(x0) # deals with arrays and scalars
    a = phs.(broadcast((a,b)->a.-b, Ref(x0), x), k, Ref(phs_exponents), Ref(diff_orders))
    p = poly.(Ref(x0), poly_exponents, Ref(diff_orders))
    return a, p
end


end # module
