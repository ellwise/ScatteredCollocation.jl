using ScatteredCollocation
using JuAFEM
using Tensors
using Plots
using NearestNeighbors
using IterativeSolvers
using SparseArrays
using LinearAlgebra

hex = false
geoshape = hex ? Hexahedron : Tetrahedron
refshape = hex ? RefCube    : RefTetrahedron
order = hex ? 2 : 1;

numdims = 3
corner1 = Vec{numdims}((0.0, 0.0, 0.0))
corner2 = Vec{numdims}((10.0, 1.0, 1.0))
grid = generate_grid(geoshape, (60, 6, 6), corner1, corner2)

p = [Array(getcoordinates(node)) for node in getnodes(grid)]
num_points = length(p)
tree = KDTree(hcat(p...); leafsize=10)
nodesets = Dict(:left=>Set{Int}(), :interior=>Set{Int}(), :right=>Set{Int}())
for (cell, face) in getfaceset(grid, "left")
    for node_idx in JuAFEM.faces(grid.cells[cell])[face]
        push!(nodesets[:left], node_idx)
    end
end
for (cell, face) in getfaceset(grid, "right")
    for node_idx in JuAFEM.faces(grid.cells[cell])[face]
        push!(nodesets[:right], node_idx)
    end
end
all_nodes = Set{Int}(collect(1:num_points))
nodesets[:interior] = setdiff(all_nodes, union(nodesets[:left], nodesets[:right]))




E = 200e9
ν = 0.3
λ = E*ν / ((1+ν) * (1 - 2ν))
μ = E / (2(1+ν))
F = 1e8 # y component
D = 1.925

# displacement formulation
# https://www.win.tue.nl/casa/meetings/seminar/previous/_abstract060524_files/seminar-24-5-6.pdf
# http://hplgit.github.io/scaling-book/doc/pub/book/html/._scaling-book-solarized010.html
# https://en.wikipedia.org/wiki/Linear_elasticity#Displacement_formulation
# Navier equation: (λ+μ).grad(div(u)) + μ.vlap(u) + ρ.b = ρ Dt.Dt.u
# Homogenenous Dirichlet: u = 0
# Traction: λ.div(u).n + μ.dot((u.div(n) + grad(u), n) = F
# Elastostatics: Dt.Dt.u = 0

# In my case n = [1, 0, 0]

# (λ+2μ).grad(div(u)) - μ.curl(curl(u)) + F = 0
# (λ+μ).grad(div(u)) + μ.vlap(u) + F = 0
# div(u) = Dx.u1 + Dy.u2 + Dz.u3
# grad(div(u)) = [Dx.Dx.u1 + Dx.Dy.u2 + Dx.Dz.u3,
#                 Dy.Dx.u1 + Dy.Dy.u2 + Dy.Dz.u3,
#                 Dz.Dx.u1 + Dz.Dy.u2 + Dz.Dz.u3]
# vlap(u) = [lap(u1), lap(u2), lap(u3)]
# lap(un) = (Dx.Dx + Dy.Dy + Dz.Dz) * un
# L = Dx.Dx + Dy.Dy + Dz.Dz
# A1 = [Dx.Dx Dx.Dy Dx.Dz;
#       Dy.Dx Dy.Dy Dy.Dz;
#       Dz.Dx Dz.Dy Dz.Dz]
# A2 = [L 0 0;
#       0 L 0;
#       0 0 L]
# A = (λ+μ).A1 + μ.A2
# Linear operator needs [2 0 0], [0 2 0], [0 0 2], [1 1 0], [1 0 1], [0 1 1]
# Boundary conditions need [0 0 0]

# SPECIFY DISPLACEMENT ON LEFT BOUNDARY, ZERO FORCE ON INTERIOR, GIVEN FORCE ON RIGHT
# RHS is force and boundary condition displacement

# in practice, I would need to define normals at each node from face normals

numvars = numdims # three displacements
L = zeros(Float64, numvars*num_points, numvars*num_points)
B = zeros(Float64, numvars*num_points, numvars*num_points)
rhs = zeros(Float64, numvars*num_points)
for row in nodesets[:left]
    p0 = p[row]
    highest_order = 0
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, p0, stencil_size)
    R = ScatteredCollocation.differentiate(p0, p[nn], [0, 0, 0])
    L[row+0*num_points, nn.+0*num_points] = R
    L[row+1*num_points, nn.+1*num_points] = R
    L[row+2*num_points, nn.+2*num_points] = R
end
for row in nodesets[:right]
    p0 = p[row]
    highest_order = 0
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, p0, stencil_size)
    R = ScatteredCollocation.differentiate(p0, p[nn], [0, 0, 0])
    L[row+0*num_points, nn.+0*num_points] = R
    L[row+1*num_points, nn.+1*num_points] = R
    L[row+2*num_points, nn.+2*num_points] = R
    rhs[row+1*num_points] = D
end
for row in nodesets[:interior]
    p0 = p[row]
    highest_order = 2
    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
    nn, _ = knn(tree, p0, stencil_size)
    DxDx = ScatteredCollocation.differentiate(p0, p[nn], [2, 0, 0])
    DxDy = ScatteredCollocation.differentiate(p0, p[nn], [1, 1, 0])
    DxDz = ScatteredCollocation.differentiate(p0, p[nn], [1, 0, 1])
    DyDx = DxDy
    DyDy = ScatteredCollocation.differentiate(p0, p[nn], [0, 2, 0])
    DyDz = ScatteredCollocation.differentiate(p0, p[nn], [0, 1, 1])
    DzDx = DxDz
    DzDy = DyDz
    DzDz = ScatteredCollocation.differentiate(p0, p[nn], [0, 0, 2])
    # A = (λ+μ).A1 + μ.A2
    L[row+0*num_points, nn.+0*num_points] = (λ+μ).*DxDx + μ.*(DxDx+DyDy+DzDz)
    L[row+0*num_points, nn.+1*num_points] = (λ+μ).*DxDy
    L[row+0*num_points, nn.+2*num_points] = (λ+μ).*DxDz
    L[row+1*num_points, nn.+0*num_points] = (λ+μ).*DyDx
    L[row+1*num_points, nn.+1*num_points] = (λ+μ).*DyDy + μ.*(DxDx+DyDy+DzDz)
    L[row+1*num_points, nn.+2*num_points] = (λ+μ).*DyDz
    L[row+2*num_points, nn.+0*num_points] = (λ+μ).*DzDx
    L[row+2*num_points, nn.+1*num_points] = (λ+μ).*DzDy
    L[row+2*num_points, nn.+2*num_points] = (λ+μ).*DzDz + μ.*(DxDx+DyDy+DzDz)
    R = ScatteredCollocation.differentiate(p0, p[nn], [0, 0, 0])
    B[row+0*num_points, nn.+0*num_points] = R
    B[row+1*num_points, nn.+1*num_points] = R
    B[row+2*num_points, nn.+2*num_points] = R
end
#for row in nodesets[:right]
#    p0 = p[row]
#    highest_order = 1
#    stencil_size = ScatteredCollocation.smallest_stencil(numdims, highest_order)
#    nn, _ = knn(tree, p0, stencil_size)
#    Dx = ScatteredCollocation.differentiate(p0, p[nn], [1, 0, 0])
#    Dy = ScatteredCollocation.differentiate(p0, p[nn], [0, 1, 0])
#    Dz = ScatteredCollocation.differentiate(p0, p[nn], [0, 0, 1])
    # λ.div(u).n + μ.dot((u.div(n) + grad(u), n) = F ??????
    # σxx = [(2μ+λ).Dx λ.Dy λ.Dz] * u;
    # σxy = [μ.Dy μ.Dx 0] * u
    # dot(σ, n) = t
    # n = [1, 0, 0]
    # dot(σ, n) = n1*[σxx σyx σzx]'
    #           + n2*[σxy σyy σzy]'
    #           + n3*[σxz σyz σzz]'
    #           = n1*[(2μ+λ).Dx      λ.Dy      λ.Dz;
    #                      μ.Dy      μ.Dx         0;
    #                      μ.Dz         0      μ.Dx]
    #           + n2*[     μ.Dy      μ.Dx         0;
    #                      λ.Dx (2μ+λ).Dy      λ.Dz;
    #                         0      μ.Dz      μ.Dy]
    #           + n3*[     μ.Dz         0      μ.Dx;
    #                         0      μ.Dz      μ.Dy;
    #                      λ.Dx      λ.Dy (2μ+λ).Dz]
#    L[row+0*num_points, nn.+0*num_points] = (2μ+λ).*Dx
#    L[row+0*num_points, nn.+1*num_points] = λ.*Dy
#    L[row+0*num_points, nn.+2*num_points] = λ.*Dz
#    L[row+1*num_points, nn.+0*num_points] = μ.*Dy
#    L[row+1*num_points, nn.+1*num_points] = μ.*Dx
    #L[row+1*num_points, nn.+2*num_points] = 0
#    L[row+2*num_points, nn.+0*num_points] = μ.*Dz
    #L[row+2*num_points, nn.+1*num_points] = 0
#    L[row+2*num_points, nn.+2*num_points] = μ.*Dx
#    rhs[row+1*num_points] = F
#end

#u = L \ rhs
#u = gmres(L, rhs; maxiter=100, tol=1e-3, verbose=true)
#ux = u[0*num_points+1:1*num_points]
#uy = u[1*num_points+1:2*num_points]
#uz = u[2*num_points+1:3*num_points]
#x = map(y->y[1], p)
#y = map(y->y[2], p)
#z = map(y->y[3], p)

#plotly()
#scatter(x+ux, y+uy, z+uz, markersize=2)

λ = LinearAlgebra.eigvals(L, B)
λ = λ[sortperm(real(λ))]

plotly()
scatter(real(λ), imag(λ), markersize=2)
