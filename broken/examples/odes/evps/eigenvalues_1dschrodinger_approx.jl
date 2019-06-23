# Example 3 from Chebfun Blocks paper

import ScatteredCollocation
using NearestNeighbors
using Plots
using SparseArrays
using LinearAlgebra

# http://www.chebfun.org/examples/ode-eig/Eigenstates.html

num_points = 128
d = 3.0
h = 0.1
V(x) = x^2
V(x) = 10 - 10*(abs(x)<1)
V(x) = abs(x)
V(x) = sqrt(abs(x)+0.1)
V(x) = 0.5*(abs(x-0.5)<0.5)
V(x) = 0.5*exp(-2*(x-0.5)^2)

stencil_size = 2*ScatteredCollocation.smallest_stencil(1, 2)

x = collect(range(-d, stop=d, length=num_points))
nodesets = Dict(:boundary=>Set([1, num_points]), :interior=>Set(2:num_points-1))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

# JULIA IS COLUMN-MAJOR! NEED TO TRANSPOSE MY LOOPS!
ii = Vector{Int}(undef, stencil_size*num_points)
jj = Vector{Int}(undef, stencil_size*num_points)
cc = Vector{Float64}(undef, stencil_size*num_points)
bb = Vector{Float64}(undef, stencil_size*num_points)
for col = nodesets[:boundary]
    index = (col-1)*stencil_size
    x0 = x[col]
    nn, _ = knn(tree, [x0], stencil_size)
    ii[index+1:index+stencil_size] .= col
    jj[index+1:index+stencil_size] = collect(1:num_points)[nn]
    cc[index+1:index+stencil_size] = ScatteredCollocation.differentiate(x0, x[nn], [0])
    bb[index+1:index+stencil_size] .= 0.0
end
for col = nodesets[:interior]
    index = (col-1)*stencil_size
    x0 = x[col]
    nn, _ = knn(tree, [x0], stencil_size)
    ii[index+1:index+stencil_size] .= col
    jj[index+1:index+stencil_size] = collect(1:num_points)[nn]
    coeffs1 = ScatteredCollocation.differentiate(x0, x[nn], [2])
    coeffs2 = ScatteredCollocation.differentiate(x0, x[nn], [0])
    cc[index+1:index+stencil_size] = -h.*coeffs1 + V(x0).*coeffs2
    bb[index+1:index+stencil_size] = ScatteredCollocation.differentiate(x0, x[nn], [0])
end
L = sparse(ii, jj, cc)
B = sparse(ii, jj, bb)

# find and sort eigenvalues
λ = LinearAlgebra.eigvals(Array(L), Array(B))
λ = λ[sortperm(real(λ))]

# print the 6 smallest (should be 1, 3, 5, 7, 9, 11)
print(real(λ[1:10]))

# solve for and plot the smallest eigenvectors
#u = nullspace.(Ref(Array(L)) .- λ[1:6].*Ref(Array(B))) .+ real(λ[1:6])
u1 = nullspace(Array(L - real(λ[1])*B)) .+ real(λ[1])
u2 = nullspace(Array(L - real(λ[2])*B)) .+ real(λ[2])
u3 = nullspace(Array(L - real(λ[3])*B)) .+ real(λ[3])
u4 = nullspace(Array(L - real(λ[4])*B)) .+ real(λ[4])
u5 = nullspace(Array(L - real(λ[5])*B)) .+ real(λ[5])
u6 = nullspace(Array(L - real(λ[6])*B)) .+ real(λ[6])
u7 = nullspace(Array(L - real(λ[7])*B)) .+ real(λ[7])
u8 = nullspace(Array(L - real(λ[8])*B)) .+ real(λ[8])
u9 = nullspace(Array(L - real(λ[9])*B)) .+ real(λ[9])

plotly()
plot(x, [u1, u2, u3, u4, u5, u6, u7, u8, u9], markersize=2)
