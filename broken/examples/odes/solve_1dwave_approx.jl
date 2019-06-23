# Example 1 from Chebfun Blocks paper

using ScatteredCollocation
using NearestNeighbors
using Plots
using SparseArrays

num_points = 256
T = Array{Float64}
stencil_size = ScatteredCollocation.smalleststencil(1, 2)*4

x = collect(range(-1.0, stop=1.0, length=num_points))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

ind = [1, num_points]
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
Au = assemble(T, x, ind, nn, identityop)
Av = zeros(Float64, length(ind), num_points)
f = [0.0, 1.0]

ind = collect(2:num_points-1)
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
Bu = assemble(T, x, ind, nn, gradientop)
Bv = assemble(T, x, ind, nn, (y0, y) -> -10*identityop(y0, y))
gu = zeros(Float64, length(ind))

ind = collect(1:num_points)
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
Cu = assemble(T, x, ind, nn, (y0, y) -> 30*exp(3*y0)*identityop(y0, y))
Cv = assemble(T, x, ind, nn, gradientop)
gv = zeros(Float64, length(ind))

L = [Au Av; Bu Bv; Cu Cv]
rhs = [f; gu; gv]

u = L \ rhs
v = u[num_points+1:end]
u = u[1:num_points]

plotly()
plot(x, [u, v], markersize=2)
