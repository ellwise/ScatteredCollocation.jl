# Example 1 from Chebfun Blocks paper

import ScatteredCollocation
using NearestNeighbors
using Plots
using LinearAlgebra

# http://www.chebfun.org/examples/ode-eig/OrrSommerfeld.html

numdims = 1
num_points = 256
alph = 1.0
Re = 2000
#Re = 5772.22

x = collect(range(-1.0, stop=1.0, length=num_points))
nodesets = Dict(:boundary=>Set([1, num_points]), :interior=>Set(2:num_points-1))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

# JULIA IS COLUMN-MAJOR! NEED TO TRANSPOSE MY LOOPS!
A = zeros(Complex{Float64}, num_points, num_points)
B = zeros(Complex{Float64}, num_points, num_points)
for row = nodesets[:boundary]
    x0 = x[row]
    highest_order = 0
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    A[row, nn] = ScatteredCollocation.differentiate(x0, x[nn], [0])
    B[row, nn] .= 0.0
end
for row = nodesets[:interior]
    x0 = x[row]
    highest_order = 6
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    D4 = ScatteredCollocation.differentiate(x0, x[nn], [4])
    D2 = ScatteredCollocation.differentiate(x0, x[nn], [2])
    R = ScatteredCollocation.differentiate(x0, x[nn], [0])
    A[row, nn] = (D4 - 2*alph^2 .* D2 + alph^4 .* R)./Re - 2im*alph .* R - 1im*alph*(1-x0^2).*(D2 - alph^2 .* R)
    B[row, nn] = D2 - R
end

λ = LinearAlgebra.eigvals(Array(A), Array(B))
λ = λ[sortperm(real(λ))]
λ = λ[end-50:end]

# solve for and plot the smallest eigenvectors
plotly()
scatter(real(λ), imag(λ), markersize=2)
