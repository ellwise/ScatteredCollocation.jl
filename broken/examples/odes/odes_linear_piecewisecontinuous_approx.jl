# http://www.chebfun.org/examples/ode-linear/DawsonIntegral.html

using ScatteredCollocation
using NearestNeighbors
using Plots
using IterativeSolvers
using SparseArrays

numdims = 1
num_points = 200
T = SparseMatrixCSC#Array{Float64}
stencil_size = ScatteredCollocation.smalleststencil(numdims, 1)

x = collect(range(0.0, stop=8.0, length=num_points))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

ind = [1]
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
A = assemble(T, x, ind, nn, identityop)
f = zeros(Float64, length(ind))

ind = collect(2:num_points)
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
op(centre, stencil) = linearop(centre, stencil, [1, 0], [1, 1])
B = assemble(T, x, ind, nn, op)
g = sin.(x[ind].^2)

#u = [A; B] \ [f; g]
u = IterativeSolvers.gmres([A; B], [f; g]; maxiter=20, tol=1e-3, verbose=true)

plotly()
plot(x, u, markersize=2)
