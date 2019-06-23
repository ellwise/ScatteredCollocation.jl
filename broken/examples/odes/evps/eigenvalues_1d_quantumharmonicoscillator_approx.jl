# Example 1 from Chebfun Blocks paper

import ScatteredCollocation
using NearestNeighbors
using Plots
using LinearAlgebra

# http://www.chebfun.org/examples/ode-eig/OpticalResponse.html

numdims = 1
num_points = 256
L = 8.0
E = 0.0

x = collect(range(-L, stop=L, length=num_points))
nodesets = Dict(:boundary=>Set([1, num_points]), :interior=>Set(2:num_points-1))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

# JULIA IS COLUMN-MAJOR! NEED TO TRANSPOSE MY LOOPS!
L = zeros(Float64, num_points, num_points)
B = zeros(Float64, num_points, num_points)
for row = nodesets[:boundary]
    x0 = x[row]
    highest_order = 0
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    L[row, nn] = ScatteredCollocation.differentiate(x0, x[nn], [0])
    B[row, nn] .= 0.0
end
for row = nodesets[:interior]
    x0 = x[row]
    highest_order = 2
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    coeffs1 = ScatteredCollocation.differentiate(x0, x[nn], [2])
    coeffs2 = ScatteredCollocation.differentiate(x0, x[nn], [0])
    L[row, nn] = -0.5.*coeffs1 + (2*x0^2 + E*x0).*coeffs2
    B[row, nn] = coeffs2
end

λ = LinearAlgebra.eigvals(Array(L), Array(B))
λ = λ[sortperm(real(λ))]

# solve for and plot the smallest eigenvectors
#u = nullspace.(Ref(Array(L)) .- λ[1:6].*Ref(Array(B))) .+ real(λ[1:6])
u1 = nullspace(Array(L - real(λ[1])*B))
u2 = nullspace(Array(L - real(λ[2])*B))
u3 = nullspace(Array(L - real(λ[3])*B))
u4 = nullspace(Array(L - real(λ[4])*B))

plotly()
plot(x, [u1, u2, u3, u4], markersize=2)
