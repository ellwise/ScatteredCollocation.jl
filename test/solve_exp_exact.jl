using ScatteredCollocation
using NearestNeighbors
using Test
using SparseArrays
using LinearAlgebra
using IterativeSolvers
#using IncompleteLU
#using AlgebraicMultigrid
# Krylov, KrylovMethods

using AbstractPlotting

# Solve:
# dudx - u = 0
#     u(0) = 1

num_points = 100
k = 3
accuracy_order = 3
T = Array{Float64}#SparseMatrixCSC{BigFloat}#

x = collect(range(0.0, stop=1.0, length=num_points))
tree = KDTree(hcat([[y] for y in x]...); leafsize=10)

bidxs = [1]
iidxs = setdiff(1:num_points,bidxs)

# assemble interior/boundary operators and RHSs
sca = SCAssembler3(k, accuracy_order, tree)
A = derivative(T, sca, 1, iidxs) - derivative(T, sca, 0, iidxs)
B = derivative(T, sca, 0, bidxs)
f = zeros(Float64, length(iidxs))
g = ones(Float64, length(bidxs))

# project the BCs onto the interior
A = hcat(A[:,iidxs], A[:,bidxs])
B = hcat(B[:,iidxs], B[:,bidxs])
Ap, Bp, fp, gp = projbc(A, B, f, g)

#print(LinearAlgebra.cond(Array([A; B])))
#print('\n')

e = eigvals(Ap)
println(maximum(real.(e)))
ymin, ymax = extrema(imag.(e))
markersize = 0.01*(ymax-ymin)
scene = AbstractPlotting.scatter(real.(e), imag.(e), markersize=markersize)
display(scene)

# solve system, extrapolate to boundary, check solution
u = Ap\fp#IterativeSolvers.gmres(Ap, fp; verbose=true)#
v = zeros(Float64, num_points)
v[iidxs] = u
v[bidxs] = gp-Bp*u

#ymin, ymax = extrema(v)
#markersize = 0.01*(ymax-ymin)
#scene = AbstractPlotting.lines(x, exp.(x))
#AbstractPlotting.scatter!(scene, x, v, markersize=markersize)
#display(scene)

#@test isapprox(u, exp.(x); rtol=sqrt(1e-3), atol=0)
#@test u ≈ exp.(x)
