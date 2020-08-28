using ScatteredCollocation
using Makie
using LinearAlgebra
using Arpack
using SparseArrays, NearestNeighbors
using DifferentialEquations, OrdinaryDiffEq, DiffEqOperators, DiffEqCallbacks
using LabelledArrays
SC = ScatteredCollocation

meshfile = joinpath(@__DIR__,"../../meshes/dolfin_fine.xml")
k = 5
accuracy_order = 3
tol = 1e-6
ε = 1e-1
Δt = 0.01
Δx = 0.05
tspan = (0.0, 2.0)
T1 = Float64
T2 = SparseMatrixCSC{T1}

# generate collocation points
#dolfin_mesh = ScatteredCollocation.readdolfinmesh(meshfile)
#x = [[y...] for y in dolfin_mesh.vertices]
#num_points = length(x)
#is_boundary = ScatteredCollocation.getboundary(dolfin_mesh)
#x, is_boundary = ScatteredCollocation.griddedSquarePoints(Δx/2); x *= 2
x, is_boundary = SC.concentricDiscPoints(Δx); x = x.+Ref([0.5, 0.5])
npts = length(x)
tree = KDTree(hcat(map(y->[y...], x)...))
idxs = collect(1:npts)
bidxs = idxs[is_boundary]
iidxs = setdiff(idxs, bidxs)

# reorder DOFs (not certin this is correct)
# doesn't diagonalise particularly well...
#permute!(x, tree.indices)
#bidxs = invperm(tree.indices)[bidxs]
#iidxs = setdiff(idxs, bidxs)
#tree = KDTree(hcat(map(y->[y...], x)...))

# assemble interior/boundary operators and define RHSs
sca = Assembler(k, accuracy_order, tree)
E  = derivative(T2, sca, [0, 0], bidxs)
∂x = derivative(T2, sca, [1, 0], idxs)
∂y = derivative(T2, sca, [0, 1], idxs)
∇² = derivative(T2, sca, [2, 0], idxs) + derivative(T2, sca, [0, 2], idxs)
∇⁴ = -derivative(T2, sca, [4, 0], idxs) - 2*derivative(T2, sca, [2, 2], idxs) - derivative(T2, sca, [0, 4], idxs)
#∇⁶ =   derivative(T2, sca, [6, 0], idxs) + 3*derivative(T2, sca, [4, 2], idxs)
#   + 3*derivative(T2, sca, [2, 4], idxs) +   derivative(T2, sca, [0, 6], idxs)
E = zero(∂x) + I
E = E[bidxs,:]
Z1 = zero(∂x)
Z2 = zero(E)
#A = [Z1 Z1 -Dx; Z1 Z1 -Dy; -Dx[iidxs,:] -Dy[iidxs,:] Z1[iidxs,:]]
A = [Z1 Z1 -∂x; Z1 Z1 -∂y; -∂x[iidxs,:] -∂y[iidxs,:] ε*∇⁴[iidxs,:]]
B = [Z2 Z2 E]
#B1  = ∂x[bidxs,:] .* map(y->y[1]-0.5, x[bidxs]) + ∂y[bidxs,:] .* map(y->y[2]-0.5, x[bidxs])
#B1 = -∂x[bidxs,:] .* [y[1]-0.5 for y in x[bidxs]] -∂y[bidxs,:] .* [y[2]-0.5 for y in x[bidxs]]
#B1 = B1 ./ [norm(y) for y in x[bidxs]]
#B = [Z2 Z2 B1]
f = zeros(T1, 2*length(idxs)+length(iidxs))
g = zeros(T1, length(bidxs))

# define initial condition
x0 = [0.75, 0.8]
p0 = exp.(-norm.(x.-Ref(x0),2).^2 ./ (0.1^2)) # pressure
vx0 = zero(p0)                    # velocity
vy0 = zero(p0)
#u0 = @LArray [vx0; vy0; p0[iidxs]] (vx=1:npts, vy=(npts+1):(2*npts), p=(2*npts+1):(2*npts+length(iidxs)))
u0 = [vx0; vy0; p0[iidxs]]

# project the BCs onto the interior
A = A[:, vcat(idxs, idxs.+npts, iidxs.+2*npts, bidxs.+2*npts)]
B = B[:, vcat(idxs, idxs.+npts, iidxs.+2*npts, bidxs.+2*npts)]
𝓐, 𝓑, 𝓯, 𝓰 = projbc(A, B, f, g)
dropzeros!(𝓐)
dropzeros!(𝓑)
print("Sparsity = $(1-nnz(𝓐)/prod(size(𝓐)))")

p = symrcm(𝓐)
ip = invperm(p)
𝓐 = 𝓐[p, p]
𝓑 = 𝓑[:,p]
u0 = u0[p]
"""
fact = svd(Array(𝓐); full=true)
fact.S[fact.S .< 1e-1*maximum(fact.S)] .= 0
𝓐 = fact.U*Diagonal(fact.S)*fact.Vt
"""

# plot eigenvalues
#λ, _ = Arpack.eigs(∇⁴; nev=npts, which=:SM)
#fact = svd(Array(𝓐); full=true)
#idxs = fact.S .> tol*maximum(fact.S)
#Ai = fact.V[:,idxs] * Diagonal(1 ./ fact.S[idxs]) * transpose(fact.U)[idxs,:]
"""
λ = eigvals(Array(𝓐))# + 10000*Ai)
scene = AbstractPlotting.Scene(scale_plot=true)
scatter!(real.(λ), imag.(λ), markersize=1000)#, xlim=(-0.015,0.015), ylim=(-2,2))
display(scene)
afdsfadsfads
"""

#scene = AbstractPlotting.Scene(scale_plot=false)
#scspy!(scene, 𝓐)
#display(scene)
#fasdafdsfsadfsd

# initialise the plot and define the update function
p = Node(p0)
clims = 0.2*maximum(abs.(p0))*[-1, 1]
sortres = true
_, ds = knn(tree, SC.getpoints(tree), 2, sortres)
Δxs = [d[2] for d in ds]
scene = AbstractPlotting.Scene(scale_plot=false)
scplot!(scene, tree, p; colorrange=clims, colormap=:pu_or)
display(scene)

# define the ODE problem and solve
function plotter(integrator)
    pnew = copy(p0)
    pnew[iidxs] = integrator.u[ip[(2*npts+1):end]]
    pnew[bidxs] = 𝓰 - 𝓑*integrator.u
    push!(p, pnew)
    sleep(1/60)
end
cb = PeriodicCallback(plotter, Δt; save_positions=(false,false))
𝓞 = DiffEqArrayOperator(𝓐)
prob = ODEProblem(𝓞, u0, tspan)
timekwargs  = Dict(:save_everystep=>false, :save_start=>false)
tolkwargs   = Dict(:reltol=>tol, :abstol=>tol)
otherkwargs = Dict(:alg_hints=>[:memorybound, :auto], :callback=>cb)
u = solve(prob; timekwargs..., tolkwargs..., otherkwargs...)

# extract components
p = copy(p0)
p[iidxs] = u[ip[(2*npts+1):end],end]
p[bidxs] = 𝓰 - 𝓑*u[:,end]
vx = u[ip[idxs],end]
vy = u[ip[idxs.+npts],end]