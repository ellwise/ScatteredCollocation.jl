# http://www.chebfun.org/examples/ode-linear/DawsonIntegral.html

import ScatteredCollocation
using NearestNeighbors
using Plots

numdims = 1
num_points = 257 # deliberately odd to ensure a point at x=0
T = Array{Float64}
stencil_size = ScatteredCollocation.smalleststencil(numdims, 1)

x = collect(range(-5.0, stop=5.0, length=num_points))
nodesets = Dict(:boundary=>Set{Int}(), :interior=>Set{Int}())
push!(nodesets[:boundary], floor(num_points/2.0)+1)
all_nodes = Set{Int}(collect(1:num_points))
nodesets[:interior] = setdiff(all_nodes, nodesets[:boundary])
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

ind = collect(nodesets[:boundary])
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
A = assemble(T, x, ind, nn, identityop)
f = zeros(Float64, length(ind))

ind = collect(nodesets[:interior])
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
op(centre, stencil) = linearop(centre, stencil, [1, 0], [1, 2*centre])
B = assemble(T, x, ind, nn, op)
g = ones(Float64, length(ind))

u = [A; B] \ [f; g]

plotly()
plot(x, u, markersize=2)
