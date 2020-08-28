using ScatteredCollocation
using NearestNeighbors
using Test
using SparseArrays
using IterativeSolvers
using AbstractPlotting
#using IncompleteLU
#using AlgebraicMultigrid
# Krylov, KrylovMethods

# Solve:      dudx - u = 0
# Subject to:     u(0) = 1
# Answer:            u = exp(x)

num_points = 1000
xspan = (0.0, 1.0)
k = 3
order = 7
T = SparseMatrixCSC{Float64}#Array{Float64}#

# generate collocation points
x = collect(range(xspan[1], stop=xspan[2], length=num_points))
tree = KDTree(hcat([[y] for y in x]...); leafsize=10)
bidxs = [1]
iidxs = setdiff(1:num_points, bidxs)

# assemble interior/boundary operators and RHSs
sca = SCAssembler3(k, order, tree)
A = derivative(T, sca, 1, iidxs) - derivative(T, sca, 0, iidxs)
B = derivative(T, sca, 0, bidxs)
f = zeros(Float64, length(iidxs))
g = ones(Float64, length(bidxs))

# project the BCs onto the interior
A = hcat(A[:,iidxs], A[:,bidxs])
B = hcat(B[:,iidxs], B[:,bidxs])
Ap, Bp, fp, gp = projbc(A, B, f, g)

# solve system, extrapolate to boundary, check solution
u0 = ones(length(iidxs)) # initial guess
u = IterativeSolvers.bicgstabl!(u0, Ap, fp; max_mv_products=10*num_points, verbose=false)
#u = Ap\fp#
v = Array{Float64,1}(undef, num_points)
v[iidxs] = u
v[bidxs] = gp-Bp*u

#ymin, ymax = extrema(v)
#markersize = 0.01*(ymax-ymin)
#scene = AbstractPlotting.lines(x, exp.(x))
#AbstractPlotting.scatter!(scene, x, v, markersize=markersize)
#display(scene)

@test v .+ 1 ≈ exp.(x) .+ 1
