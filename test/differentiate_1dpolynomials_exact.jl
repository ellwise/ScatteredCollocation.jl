import ScatteredCollocation
using Test
using IterTools

"""
One-dimensional differentiation of monomials, structured nodes, multiple orders
of accuracy, multiple orders of derivative.
"""

numdims = 1
num_eval = 100
highest_order = 1

# differentiate on and between nodes
eval_points = range(-2, stop=2, length=num_eval)

min_num_nodes = ScatteredCollocation.smallest_stencil(numdims, highest_order)

# test varying stencil sizes
for num_nodes = min_num_nodes:2:round(Int, num_eval/10)
    colloc_nodes = collect(range(-2, stop=2, length=num_nodes))

    # should be exact for all polynomials up to a given order
    num_poly = floor(Int, num_nodes/2)
    #maxorder = prod(IterTools.nth(ScatteredCollocation.monomial_exponents(), num_poly))
    for monomial_index = 1:num_poly

        # generate samples and true derivative
        exponents = IterTools.nth(ScatteredCollocation.monomial_exponents(numdims), monomial_index)
        u = ScatteredCollocation.poly.(colloc_nodes, exponents)
        diff_order = 1
        up_true = ScatteredCollocation.poly.(eval_points, Ref(exponents), Ref(diff_order))

        # assemble differentiation matrix
        D = zeros(Float64, num_eval, num_nodes)
        for row = 1:num_eval
            x0 = eval_points[row]
            D[row, :] = ScatteredCollocation.differentiate(x0, colloc_nodes, diff_order)
        end

        # differentiate samples and compare
        up_approx = D*u
        @test up_approx .+ 1 ≈ up_true .+ 1 # 0.0 ≈ 0.0 returns false!

    end # poly_order
end # num_nodes
