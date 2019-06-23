
import ScatteredCollocation
using NearestNeighbors
using Plots
using IterativeSolvers
using SparseArrays
using LinearAlgebra
using JuAFEM

numdims = 2

f = 1
k = 1

grid = JuAFEM.generate_grid(Quadrilateral, (20,20))
p = [Array(getcoordinates(node)) for node in getnodes(grid)]
num_points = length(p)
tree = KDTree(hcat(p...); leafsize=10)
facesets = union(getfaceset(grid, "left"), getfaceset(grid, "right"), getfaceset(grid, "top"), getfaceset(grid, "bottom"))
nodesets = Dict(:boundary=>Set{Int}(), :interior=>Set{Int}())
for (cell, face) in facesets
    for node_idx in JuAFEM.faces(grid.cells[cell])[face]
        push!(nodesets[:boundary], node_idx)
    end
end
all_nodes = Set{Int}(collect(1:num_points))
nodesets[:interior] = setdiff(all_nodes, nodesets[:boundary])

L = zeros(Float64, num_points, num_points)
rhs = zeros(Float64, num_points)
for row in nodesets[:boundary]
    p0 = p[row]
    highest_order = 0 # FAILS IF THIS IS ZERO!
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, p0, stencil_size)
    L[row, nn] = ScatteredCollocation.differentiate(p0, p[nn], [0, 0])
end
for row in nodesets[:interior]
    p0 = p[row]
    highest_order = 2
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, p0, stencil_size)
    DxDx = ScatteredCollocation.differentiate(p0, p[nn], [2, 0])
    DyDy = ScatteredCollocation.differentiate(p0, p[nn], [0, 2])
    L[row, nn] = DxDx + DyDy
    rhs[row] = 1
end

# below can be slower, but improves conditioning? or regularises?
u = IterativeSolvers.gmres(L, rhs; maxiter=1000, tol=1e-3, verbose=true)
x = map(y->y[1], p)
y = map(y->y[2], p)

#u = IterativeSolvers.bicgstabl(L, rhs; max_mv_products=100, tol=1e-9, verbose=true)

plotly()
scatter(x, y, u, markersize=2)

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
