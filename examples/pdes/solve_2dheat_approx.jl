
import ScatteredCollocation
using NearestNeighbors
using Plots
using IterativeSolvers
using SparseArrays
using LinearAlgebra

point_spacing = 0.03
accuracy_order = 4 # this is roughly the limit before it goes haywire
num_poly = ScatteredCollocation.unique_monomials(2, accuracy_order)
stencil_size = 2*num_poly

f = 1
k = 1

x, bnd = ScatteredCollocation.griddedSquarePoints(point_spacing)
num_points = length(x)
tree = KDTree(hcat(x...); leafsize=10)

# JULIA IS COLUMN-MAJOR! NEED TO TRANSPOSE MY LOOPS!
ii = Vector{Int}(undef, stencil_size*num_points)
jj = Vector{Int}(undef, stencil_size*num_points)
dx = Vector{Float64}(undef, stencil_size*num_points)
dy = Vector{Float64}(undef, stencil_size*num_points)

for col = 1:num_points
    index = (col-1)*stencil_size
    x0 = x[col]
    nn, _ = knn(tree, x0, stencil_size)
    ii[index+1:index+stencil_size] .= col
    jj[index+1:index+stencil_size] = collect(1:num_points)[nn]
    dx[index+1:index+stencil_size] = ScatteredCollocation.differentiate(x0, x[nn], [1, 0])
    dy[index+1:index+stencil_size] = ScatteredCollocation.differentiate(x0, x[nn], [0, 1])
end
Dx = sparse(ii, jj, dx)
Dy = sparse(ii, jj, dy)
L = - (Dx*(k .* Dx) + Dy*(k .* Dy))

# add boundary condition
ii, jj, _ = findnz(L)
vv = nonzeros(L)
for m = 1:length(vv)
    if bnd[ii[m]]
        if bnd[jj[m]]
            vv[m]  = 1.0
        else
            vv[m] = 0.0
        end
    end
end

# add RHS and solve
rhs = f .* ones(Float64, num_points)
rhs[bnd] .= 0.0
#u = L \ rhs
# below can be slower, but improves conditioning? or regularises?
u = IterativeSolvers.gmres(L, rhs; maxiter=1000, tol=1e-3, verbose=true)
#u = IterativeSolvers.bicgstabl(L, rhs; max_mv_products=100, tol=1e-9, verbose=true)

plotly()
scatter(map(y->y[1], x), map(y->y[2], x), u, markersize=2)

#λ_max, _ = powm(L)

#λ = LinearAlgebra.eigvals(Array(L))

#print(λ_max)
#print('\n')
#print(minimum(real.(λ)))
#print('\n')

#plotly()
#scatter(real.(λ), imag.(λ), markersize=2)

#print(maximum(abs.(u)))
#print('\n')
