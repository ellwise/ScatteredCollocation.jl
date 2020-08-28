using ScatteredCollocation
using NearestNeighbors
using Test
using LinearAlgebra

# Chebyshev differential equations
# https://en.wikipedia.org/wiki/Chebyshev_polynomials
# Boundary conditions:
# https://en.wikipedia.org/wiki/Chebyshev_polynomials#Roots_and_extrema

using Plots

num_points1 = 1000
num_points2 = round(Int, 0.99*num_points)
stencil_size = 20
n = 40
kind = 1

x1 = collect(range(-1.0, stop=1.0, length=num_points1))
x2 = collect(range(-1.0, stop=1.0, length=num_points2))
tree1 = KDTree(hcat([[y] for y=x1]...); leafsize=10)
tree2 = KDTree(hcat([[y] for y=x2]...); leafsize=10)

A = zeros(Float64, num_points2, num_points1)
B = zeros(Float64, 1, num_points1)
rhsA = zeros(Float64, num_points2, 1)
rhsB = zeros(Float64, 1, 1)

nn, _ = knn(tree1, reshape(x2, 1, num_points2), stencil_size)
for row = 1:num_points2
    id, ddx, d2dx2 = differentiate(x2[row], x1[nn[row]], [0, 1, 2])
    if kind==1
        A[row, nn[row]] = (1-x2[row]^2)*d2dx2 - x2[row]*ddx + (n^2)*id
    elseif kind==2
        A[row, nn[row]] = (1-x2[row]^2)*d2dx2 - 3*x2[row]*ddx + n*(n+2)*id
    end
end
nn, _ = knn(tree1, reshape(x1, 1, num_points1), stencil_size)
B[end, nn[end]] = identityop(x1[end], x1[nn[end]])
if kind==1
    rhsB[end] = 1.0
elseif kind==2
    rhsB[end] = n+1.0
end

L = [2*A'*A B'; B zeros(Float64, 1, 1)]
rhs = [2*A'*rhsA; rhsB]

u = L\rhs
u = u[1:num_points]

if kind==1
    u_true = cos.(n*acos.(x))
elseif kind==2
    u_true = sin.((n+1)*acos.(x))./sin.(acos.(x))
    # end points don't evaluate nicely
    u_true[1] = (n+1.0)*(-1)^n
    u_true[end] = n+1.0
end

plotly()
plot(x, u)
plot!(x, u_true, markerstrokestyle=:dash)

#@test isapprox(u, u_true; rtol=1e-3, atol=1e-3)
