# Example 1 from Chebfun Blocks paper

import ScatteredCollocation
using NearestNeighbors
using Plots
using LinearAlgebra

# http://www.chebfun.org/examples/ode-eig/WaveDecay.html

numdims = 1
num_points = 256
L = pi/2
a = 0.2
middle(x) = abs(x) ≤ a

x = collect(range(-L, stop=L, length=num_points))
nodesets = Dict(:boundary=>Set([1, num_points]), :interior=>Set(2:num_points-1))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

# JULIA IS COLUMN-MAJOR! NEED TO TRANSPOSE MY LOOPS!
L = zeros(Float64, num_points, num_points)
B = zeros(Float64, num_points, num_points)
rhs = Vector{Float64}(undef, num_points)
for row = nodesets[:boundary]
    x0 = x[row]
    highest_order = 0
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    L[row, nn] = ScatteredCollocation.differentiate(x0, x[nn], [0])
    B[row, nn] .= 0.0
    rhs[row] = 0.0
end
for row = nodesets[:interior]
    x0 = x[row]
    highest_order = 2
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, [x0], stencil_size)
    coeffs1 = ScatteredCollocation.differentiate(x0, x[nn], [2])
    coeffs2 = ScatteredCollocation.differentiate(x0, x[nn], [1])
    coeffs3 = ScatteredCollocation.differentiate(x0, x[nn], [0])
    L[row, nn] = coeffs1# + (2/a)*middle(x0).*coeffs2
    B[row, nn] = coeffs3
    rhs[row] = 0.0
end

λ = LinearAlgebra.eigvals(Array(L), Array(B))
λ = λ[sortperm(real(λ))]
λ = λ[end:-1:1]

print(λ[[1, 2, 20, 40]])

# solve for and plot the smallest eigenvectors
#u = nullspace.(Ref(Array(L)) .- λ[1:6].*Ref(Array(B))) .+ real(λ[1:6])
u1 = nullspace(Array(L - real(λ[1])*B))
u2 = nullspace(Array(L - real(λ[2])*B))
u20 = nullspace(Array(L - real(λ[20])*B))
u40 = nullspace(Array(L - real(λ[40])*B))

u1 = u1 ./ maximum(u1)
u2 = u2 ./ maximum(u2)
u20 = u20 ./ maximum(u20)
u40 = u40 ./ maximum(u40)

plotly()
plot(x, [u1, u2, u20, u40], markersize=2)
#scatter(real(λ), imag(λ))
