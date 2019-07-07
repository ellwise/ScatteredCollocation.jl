using ScatteredCollocation
using NearestNeighbors
using Test
using LinearAlgebra

# Chebyshev differential equations
# https://en.wikipedia.org/wiki/Chebyshev_polynomials
# Boundary conditions:
# https://en.wikipedia.org/wiki/Chebyshev_polynomials#Roots_and_extrema

using Plots

num_points = 1000
stencil_size = 20
n = 50
kind = 1

x = collect(range(-1.0, stop=1.0, length=num_points))
tree = KDTree(hcat([[y] for y=x]...); leafsize=10)

L = zeros(Float64, num_points, num_points)
rhs = zeros(Float64, num_points)

nn, _ = knn(tree, reshape(x, 1, num_points), stencil_size)
for row = 1:num_points-1
    id, ddx, d2dx2 = differentiate(x[row], x[nn[row]], [0, 1, 2])
    if kind==1
        L[row, nn[row]] = (1-x[row]^2)*d2dx2 - x[row]*ddx + (n^2)*id
    elseif kind==2
        L[row, nn[row]] = (1-x[row]^2)*d2dx2 - 3*x[row]*ddx + n*(n+2)*id
    end
end
L[end, nn[end]] = identityop(x[end], x[nn[end]])
if kind==1
    rhs[end] = 1.0
elseif kind==2
    rhs[end] = n+1.0
end

u = L\rhs

if kind==1
    u_true = cos.(n*acos.(x))
elseif kind==2
    u_true = sin.((n+1)*acos.(x))./sin.(acos.(x))
    # end points don't evaluate nicely
    u_true[1] = (n+1.0)*(-1)^n
    u_true[end] = n+1.0
end

#plotly()
#plot(x, u)
#plot!(x, u_true, markerstrokestyle=:dash)

@test isapprox(u, u_true; rtol=1e-3, atol=1e-3)
