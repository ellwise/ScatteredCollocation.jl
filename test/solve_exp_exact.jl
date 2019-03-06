using ScatteredCollocation
using NearestNeighbors
using Test
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid

# Krylov, KrylovMethods

# Solve:
# dudx - u = 0
#     u(0) = 1

num_points = 1400
stencil_size = 20
T = SparseMatrixCSC{BigFloat}#Array{Float64}#

x = collect(range(0.0, stop=1.0, length=num_points))
tree = KDTree(hcat([[y] for y in x]...); leafsize=10)

# assemble boundary condition operator
ind = [1];
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
A = assemble(T, x, ind, nn, identityop)
f = ones(Float64, length(ind))

# assemble main differential operator
op(centre, stencil) = linearop(centre, stencil, [1, 0], [1, -1])
ind = collect(2:num_points)
nn, _ = knn(tree, reshape(x[ind], 1, length(ind)), stencil_size)
B = assemble(T, x, ind, nn, op)
g = zeros(Float64, length(ind))

#print(LinearAlgebra.cond(Array([A; B])))
#print('\n')

# complete and solve system, check solution
#u = [A; B] \ [f; g]
#p = IncompleteLU.ilu([A; B], τ=0.1)
#u = IterativeSolvers.gmres([A; B], [f; g]; maxiter=20, tol=1e-12, verbose=true, Pl=p)
ml = ruge_stuben([A; B])
u = solve(ml, [f; g])
@test isapprox(u, exp.(x); rtol=sqrt(1e-3), atol=0)
#@test u ≈ exp.(x)
