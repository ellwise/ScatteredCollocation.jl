using ScatteredCollocation
using Makie
using LinearAlgebra
using Arpack
using SparseArrays, NearestNeighbors
using DifferentialEquations, OrdinaryDiffEq, DiffEqOperators, DiffEqCallbacks
using LabelledArrays
using Printf
SC = ScatteredCollocation

accuracy_order = 3
Δt = 0.03
Δx = 0.03
tspan = (0.0, 2.0)
T1 = Float64
solver = DPRKN6()
timekwargs  = Dict(:save_everystep=>false, :save_start=>false)
tolkwargs   = Dict(:reltol=>1e-6, :abstol=>1e-6)
otherkwargs = Dict(:alg_hints=>[:memorybound, :auto])

# generate collocation points
x, is_boundary = SC.concentricDiscPoints(Δx); x = x.+Ref([0.5, 0.5])
normals = [y-[0.5, 0.5] for y in x]
npts = length(x)
tree = KDTree(reduce(hcat,x))
idxs = collect(1:npts)
bidxs = idxs[is_boundary]
iidxs = setdiff(idxs, bidxs)

# assemble differential operator matrices
∂²x, ∂²y = build(T1, tree, iidxs, [[2,0],[0,2]], accuracy_order)
∂x, ∂y   = build(T1, tree, bidxs, [[1,0],[0,1]], accuracy_order)
id,      = build(T1, tree, bidxs, [[0,0]],       accuracy_order)
∇² = ∂²x + ∂²y
soft = - ∂x.*[n[1] for n in normals[bidxs]] - ∂y.*[n[2] for n in normals[bidxs]]
hard = id

# project the boundary condition onto the interior
# - define the interior pde and boundary extrapolation functions
#
# CHANGE PROJBC TO USE AFFINE OPERATORS, AND TO SEPARATELY TAKE A, B, bidxs
# ALSO, HAVE IT RETURN A WAY OF MODIFYING, RATHER THAN THE MODIFIED (I.E. THE GAUSS-JORDAN MATRIX, OR SIMILAR)
#
A = AffineOperator(∇²)
B = AffineOperator(hard)
C, D = projbc(A, B, bidxs)
E = AffineOperator(spzeros(npts,npts))
E.A[iidxs,iidxs] = C.A
E.A[bidxs,iidxs] = D.A*C.A

println("Sparsity = $(@sprintf("%d",100*(1-fnz(E.A))))%")
pde(u) = E(u)


# plot eigenvalues
λ = eigvals(Array(E.A))
scene = AbstractPlotting.Scene(scale_plot=true)
scatter!(real.(λ), imag.(λ), markersize=100)#, xlim=(-0.015,0.015), ylim=(-2,2))
display(scene)
afdsfadsfads


# define initial condition
x0 = [0.75, 0.8]
u0 = exp.(-norm.(x.-Ref(x0),2).^2 ./ (0.1^2)) # pressure

# initialise the plot and define the update function
p = Node(u0)
clims = 0.2*maximum(abs.(u0))*[-1, 1]
scene = AbstractPlotting.Scene(scale_plot=false)
scplot!(scene, tree, p; colorrange=clims, colormap=:pu_or)
display(scene)

# define the ODE problem and solve
plotter(integrator) = begin
    pnew = integrator.u[2,:]
    push!(p, integrator.u[2,:])
    sleep(1/60)
end
cb = PeriodicCallback(plotter, Δt; save_positions=(false,false))
prob = SecondOrderODEProblem((du, u, p, t) -> pde(u), zero(u0), u0, tspan)
u = solve(prob, solver; callback=cb, timekwargs..., tolkwargs..., otherkwargs...)
#p[iidxs] = u[:]
#p[bidxs] = bnd(u[:])