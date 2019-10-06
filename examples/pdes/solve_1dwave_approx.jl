using ScatteredCollocation, NearestNeighbors, SparseArrays
using DifferentialEquations, OrdinaryDiffEq, DiffEqOperators, DiffEqCallbacks
using Makie
using Sundials
using LinearAlgebra
using Printf
#using RecursiveArrayTools

npts = 200
xspan = (-2, 2)
tspan = (0.0, 4.0)
Δt = 0.05 # for plotting
tol = 1e-6
accuracy_order = 3
T1 = Float64

# generate collocation points
x = collect(range(xspan[1], stop=xspan[2], length=npts))
tree = KDTree(reduce(hcat,x))
lbidx = [1]
rbidx = [npts]
iidxs = setdiff(1:npts, [lbidx; rbidx])

# assemble interior/boundary operators and define RHSs
∇², = build(T1, tree, iidxs, [[2]], accuracy_order)
∂ₓ, = build(T1, tree, rbidx, [[1]], accuracy_order)
id, = build(T1, tree, lbidx, [[0]], accuracy_order)

# assemble the system matrix and RHS
L = spzeros(T1, npts, npts)
L[iidxs,:] = ∇²
L[lbidx,:] = id
L[rbidx,:] = ∂ₓ
y = zeros(T1, npts)

# project the boundary condition onto the interior
# - define the interior pde and boundary extrapolation functions
A, B, f, g = projbc(L, y, [lbidx; rbidx])
println("Sparsity = $(@sprintf("%d",100*(1-fnz(A))))%")
pde(u) =  A*u - f
bnd(u) = -B*u + g

# define initial condition
p0 = exp.(-(x.+1).^2 ./ (0.2^2)) # pressure
u0 = p0[iidxs]

# initialise the plot and define the update function
p = Node(p0)
scene = AbstractPlotting.Scene()
scplot!(scene, tree, p; color=:blue)
display(scene)
function plotter(integrator)
    pnew = copy(p0)
    pnew[iidxs] = integrator.u[2,:]
    pnew[[lbidx; rbidx]] = bnd(integrator.u[2,:])
    push!(p, pnew)
    sleep(1/60)
end
cb = PeriodicCallback(plotter, Δt; save_positions=(false,false))

# define the ODE problem and solve
prob = SecondOrderODEProblem((du, u, p, t) -> pde(u), zero(u0), u0, tspan)
solver = DPRKN6()
timekwargs  = Dict(:save_everystep=>false, :save_start=>false)
tolkwargs   = Dict(:reltol=>tol, :abstol=>tol)
otherkwargs = Dict(:alg_hints=>[:memorybound, :auto], :callback=>cb)
u = solve(prob, solver; timekwargs..., tolkwargs..., otherkwargs...)