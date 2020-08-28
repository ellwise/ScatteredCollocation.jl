using NearestNeighbors
using SparseArrays
using ScatteredCollocation
using IterativeSolvers
using Makie
using Arpack, LinearAlgebra
using IncompleteLU

T = Float64
k = 3
accuracy_order = 3

# parameters
P = -1000 # loading, [N]
E = 3e7   # Young's modulus [N/m2]
v = 0.3   # Poisson's ratio
D = 12    # depth of the beam [m]
L = 48    # length of the beam [m]

# dependent parameters
λ = E*v/(1+v)/(1-2v)
μ = E/2/(1+v)
λh = 2λ*μ/(λ+2μ)

# define collocation points
Δx = 0.5 # [m]
Δy = 0.5 # [m]
x = [[x,y] for x=0:Δx:L, y=-D/2:Δx:D/2]
x = x[:]
npts = length(x)
tree = KDTree(hcat(map(y->[y...], x)...))
idxs = collect(1:npts)
sidxs = [j for j=1:npts if x[j][1]==L] # static indices
fidxs = [j for j=1:npts if x[j][1]==0] # forced indices
iidxs = setdiff(idxs, vcat(sidxs, fidxs))    # interior indices

# define differential operators
∂²xx, ∂²xy, ∂²yx, ∂²yy = build(T, tree, idxs, [[2,0],[1,1],[1,1],[0,2]], accuracy_order)
∂x, ∂y = build(T, tree, idxs, [[1,0],[0,1]], accuracy_order)

# homogeneous, isotropic, linear elasticity
∇² = ∂²xx .+ ∂²yy
cauchy_navier = (λh+μ)*[∂²xx ∂²xy; ∂²yx ∂²yy] .+ μ*blockdiag(∇²,∇²)

# strain tensor
εxx = ∂²xx
εxy = 0.5*(∂²xy+∂²yx)
εyx = 0.5*(∂²yx+∂²xy)
εyy = ∂²yy

# stress tensor
σxx = λh*(εxx+εyy) + 2*μ*εxx
σxy = 2*μ*εxy
σyx = 2*μ*εyx
σyy = λh*(εxx+εyy) + 2*μ*εyy

# traction bcs on top/bottom: σ_... prescribed
# displacement bcs on left/right: u_... prescribed

# define model equations
Z = zero(E)
PDE = (λh+μ)*[∂²x ∂²xy; ∂²xy ∂²y] .+ [∇² Z; Z ∇²]
A = PDE[[iidxs; fidxs; [iidxs; fidxs].+npts], :]
B1 = [E Z][sidxs,:]
B2 = [Z E][sidxs,:]
B = [B1; B2]
f = zeros(Float64, length([iidxs; fidxs; iidxs; fidxs]))
f[length([iidxs; fidxs; iidxs])+1:end] .= P
g = zeros(Float64, 2*length(sidxs))

# project the BCs onto the interior
A = A[:, [[iidxs; fidxs]; [iidxs; fidxs].+npts; sidxs; sidxs.+npts]]
B = B[:, [[iidxs; fidxs]; [iidxs; fidxs].+npts; sidxs; sidxs.+npts]]

# plot sparsity pattern
#scene = AbstractPlotting.Scene(scale_plot=false)
#scspy!(scene, A)
#display(scene)
#fasdafdsfsadfsd

𝓐, 𝓑, 𝓯, 𝓰 = projbc(A, B, f, g)
dropzeros!(𝓐)
dropzeros!(𝓑)
println("Sparsity = $(1-nnz(𝓐)/prod(size(𝓐)))")

# reorder dofs for minimal bandwidth
p = symrcm(𝓐)
ip = invperm(p)
𝓐 = 𝓐[p,p]
𝓑 = 𝓑[:,p]
𝓯 = 𝓯[p]

# compute displacements
"""
From RBF paper II:
Unless specified otherwise in the text, the linear system of equations is solved
using the BiCGSTAB iterative method with incomplete LU factorization as preconditioner
(with 0 level of fill in) and reverse Cuthill–McKee ordering. The residual tolerance
is set to 1e-11 and the maximum number of iterations is equal to 100.
"""
LU = ilu(𝓐, τ=0.01)
#@show nnz(LU) / nnz(𝓐)
u = bicgstabl(𝓐, 𝓯; Pl=LU, tol=1e-11)#, max_mv_products=100)
u = u[ip]
𝓐 = 𝓐[p,ip]
𝓑 = 𝓑[:,ip]
𝓯 = 𝓯[ip]
ui = u
ub = 𝓰 - 𝓑*ui
vx = Array{Float64,1}(undef,npts)
vx[[iidxs; fidxs]] = ui[1:length([iidxs; fidxs])]
vx[sidxs] = ub[1:length(sidxs)]
vy = Array{Float64,1}(undef,npts)
vy[[iidxs; fidxs]] = ui[length([iidxs; fidxs])+1:end]
vy[sidxs] = ub[length(sidxs)+1:end]
v = sqrt.(vx.^2 .+ vy.^2)

# plot resulting field
scene1 = AbstractPlotting.Scene(scale_plot=false)
scplot!(scene1, tree, vx)
scene2 = AbstractPlotting.Scene(scale_plot=false)
scplot!(scene2, tree, vy)

# plot eigenvalues
#λ, _ = Arpack.eigs(𝓐; nev=npts, which=:SM)
#λ = eigvals(Array(𝓐))
#scene3 = AbstractPlotting.Scene(scale_plot=true)
#scatter!(scene3, real.(λ), imag.(λ), markersize=1000000)

scene = hbox(scene1, scene2)#, scene3)

display(scene)