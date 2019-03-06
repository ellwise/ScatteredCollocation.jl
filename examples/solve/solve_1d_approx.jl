# Example 1 from Chebfun Blocks paper

using ScatteredCollocation
using NearestNeighbors
using Plots
using SparseArrays
using IterativeSolvers

num_points = 1000
T = SparseMatrixCSC#Array{Float64}#
stencil_size = ScatteredCollocation.smalleststencil(1, 2)

x = collect(range(-1.0, stop=1.0, length=num_points))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

ind = [1, num_points]
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
A1 = assemble(SparseMatrixCSC, x, ind, nn, identityop)
A2 = assemble(Array{Float64}, x, ind, nn, identityop)
f = [0.0; 1.0]

ind = collect(2:num_points-1)
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
op(centre, stencil) = linearop(centre, stencil, [2, 0], [0.0001, centre])
B1 = assemble(SparseMatrixCSC, x, ind, nn, op)
B2 = assemble(Array{Float64}, x, ind, nn, op)

g = exp.(x[ind])

A = A1
B = B1

u = [A; B] \ [f; g]

plotly()
plot(x, u, markersize=2)
