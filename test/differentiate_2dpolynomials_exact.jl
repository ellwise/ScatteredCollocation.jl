import ScatteredCollocation
using Test
using IterTools

"""
One-dimensional differentiation of monomials, structured nodes, multiple orders
of accuracy, multiple orders of derivative.
"""

# consider using random points for evaluation, rather than upsampled?
# random isn't a good idea with tests though...
# NOTE: random points won't neccesarily include a good boundary discretisation, so best for evaluation only

numdims = 2
eval_spacing = 0.125
disc_points = ScatteredCollocation.randomDiscPoints

# differentiate on and between nodes
eval_points, _ = disc_points(eval_spacing)
num_eval = length(eval_points)

# test varying stencil sizes
for colloc_spacing = [0.5, 0.25]
    colloc_nodes, _ = disc_points(colloc_spacing)
    num_nodes = length(colloc_nodes)

    # should be exact for all polynomials up to a given order
    num_poly = floor(Int, num_nodes/2)
    for monomial_index = 1:num_poly

        # generate samples and true derivative
        exponents = IterTools.nth(ScatteredCollocation.monomial_exponents(numdims), monomial_index)
        u = ScatteredCollocation.poly.(colloc_nodes, Ref(exponents))
        ux_true = ScatteredCollocation.poly.(eval_points, Ref(exponents), Ref([1, 0]))
        uy_true = ScatteredCollocation.poly.(eval_points, Ref(exponents), Ref([0, 1]))

        # assemble differentiation matrix
        Dx = zeros(Float64, num_eval, num_nodes)
        Dy = zeros(Float64, num_eval, num_nodes)
        for row = 1:num_eval
            x0 = eval_points[row]
            Dx[row, :] = ScatteredCollocation.differentiate(x0, colloc_nodes, [1, 0])
            Dy[row, :] = ScatteredCollocation.differentiate(x0, colloc_nodes, [0, 1])
        end

        # differentiate samples and compare
        ux_approx = Dx*u
        uy_approx = Dy*u
        @test ux_approx .+ 1 ≈ ux_true .+ 1 # 0.0 ≈ 0.0 returns false!
        @test uy_approx .+ 1 ≈ uy_true .+ 1

    end # poly_order
end # num_nodes
