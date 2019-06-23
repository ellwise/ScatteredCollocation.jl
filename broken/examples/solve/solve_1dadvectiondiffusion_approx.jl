# Example 1 from Chebfun Blocks paper

using ScatteredCollocation
using NearestNeighbors
using Plots

numdims = 1
num_points = 256

x = collect(range(-1.0, stop=1.0, length=num_points))
nodesets = Dict(:left=>Set(1), :interior=>Set(2:num_points-1), :right=>Set(num_points))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

# JULIA IS COLUMN-MAJOR! NEED TO TRANSPOSE MY LOOPS!
L = zeros(Float64, num_points, num_points)
rhs = Vector{Float64}(undef, num_points)
for row = nodesets[:left]
    x0 = x[row]
    highest_order = 0
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    L[row, nn] = identityop(x0, x[nn])
    rhs[row] = 2.0
end
for row = nodesets[:right]
    x0 = x[row]
    highest_order = 1
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    L[row, nn] = differentiate(x0, x[nn], 1)
    rhs[row] = 0.0
end
for row = nodesets[:interior]
    x0 = x[row]
    highest_order = 2
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    coeffs = differentiate(x0, x[nn], [0, 1, 2])
    L[row, nn] = 0.2*coeffs[3] + coeffs[2] + coeffs[1]
    rhs[row] = exp(x0)
end

u = L \ rhs

plotly()
plot(x, u, markersize=2)
