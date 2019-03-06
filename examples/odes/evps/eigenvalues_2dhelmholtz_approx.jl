
# http://www.chebfun.org/examples/disk/Eigenfunctions.html

import ScatteredCollocation
using NearestNeighbors
using Plots
using LinearAlgebra

point_spacing = 0.1
numdims = 2
highest_order = 3
# NOTE: random points won't neccesarily include a good boundary discretisation

stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
x, bnd = ScatteredCollocation.concentricDiscPoints(point_spacing)
nodesets = Dict(:boundary=>Set(findall(bnd)), :interior=>Set(findall(.!bnd)))
num_points = length(x)
tree = KDTree(hcat(x...); leafsize=10)

L = zeros(Float64, num_points, num_points)
B = zeros(Float64, num_points, num_points)
for col = nodesets[:boundary]
    x0 = x[col]
    nn, _ = knn(tree, x0, stencil_size)
    L[col, nn] = ScatteredCollocation.differentiate(x0, x[nn], [0, 0])
    B[col, nn] .= 0.0
end
for col = nodesets[:interior]
    x0 = x[col]
    nn, _ = knn(tree, x0, stencil_size)
    d2x = ScatteredCollocation.differentiate(x0, x[nn], [2, 0])
    d2y = ScatteredCollocation.differentiate(x0, x[nn], [0, 2])
    eye = ScatteredCollocation.differentiate(x0, x[nn], [0, 0])
    L[col, nn] = d2x .+ d2y
    B[col, nn] = eye
end

λ = LinearAlgebra.eigvals(L, B)
λ = λ[sortperm(real(λ))]

print(minimum(abs.(real(λ))))

# solve for and plot the smallest eigenvectors
#u = nullspace.(Ref(Array(L)) .- λ[1:6].*Ref(Array(B))) .+ real(λ[1:6])
u = nullspace(L - minimum(abs.(real(λ)))*B)

plotly()
scatter(map(y->y[1], x), map(y->y[2], x), u, markersize=2)
#scatter(real(λ), imag(λ))
