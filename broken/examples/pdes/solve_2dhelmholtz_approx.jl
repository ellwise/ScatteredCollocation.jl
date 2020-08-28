
using ScatteredCollocation
using NearestNeighbors
using Plots
using IterativeSolvers

point_spacing = 0.1
numdims = 2
highest_order = 2
# NOTE: random points won't neccesarily include a good boundary discretisation

stencil_size = ScatteredCollocation.smalleststencil(numdims, highest_order)
x, bnd = ScatteredCollocation.concentricDiscPoints(point_spacing)
x = [(y[1], y[2]) for y in x]
nodesets = Dict(:boundary=>Set(findall(bnd)), :interior=>Set(findall(.!bnd)))
num_points = length(x)
tree = KDTree(hcat(map(y->[y...], x)...); leafsize=10)

LT = zeros(Float64, num_points, num_points)
rhs = Vector{Float64}(undef, num_points)
for col = nodesets[:boundary]
    x0 = x[col]
    nn, _ = knn(tree, [x0...], stencil_size)
    LT[nn, col] = identityop(x0, x[nn])
    rhs[col] = cos(3*atan(x0[2],x0[1]))
end
for col = nodesets[:interior]
    x0 = x[col]
    nn, _ = knn(tree, [x0...], stencil_size)
    LT[nn, col] = laplacianop(x0, x[nn])
    rhs[col] = 0.0
end
L = LT'

u = IterativeSolvers.gmres(L, rhs; maxiter=20, tol=1e-3, verbose=true)

plotly()
scatter(map(y->y[1], x), map(y->y[2], x), u, markersize=2)
