using ScatteredCollocation, NearestNeighbors, SparseArrays
using OrdinaryDiffEq, DiffEqOperators, DiffEqCallbacks
using Makie

num_points = 200
xspan = (-2, 2)
tspan = (0.0, 2.0)
Δt = 0.05 # for plotting
k = 3
accuracy_order = 3
T1 = Float32
T2 = SparseMatrixCSC{T1}#Array{T1}#

# generate collocation points
x = collect(range(xspan[1], stop=xspan[2], length=num_points))
tree = KDTree(hcat(map(y->[y...], x)...))
bidxs = [1, num_points]
iidxs = setdiff(1:num_points, bidxs)

# assemble interior/boundary operators and define RHSs
sca = SCAssembler3(k, accuracy_order, tree)
A = 0.1*derivative(T2, sca, 2, iidxs)
B =     derivative(T2, sca, 1, bidxs)
f = zeros(T1, length(iidxs))
g = zeros(T1, length(bidxs))

# define initial condition
v0 = exp.(-(x.+1).^2 ./ (0.2^2))
u0 = v0[iidxs]

# project the BCs onto the interior
A = [A[:,iidxs] A[:,bidxs]]
B = [B[:,iidxs] B[:,bidxs]]
Ap, Bp, fp, gp = projbc(A, B, f, g)
print("Sparsity = "); println(1-nnz(Ap)/prod(size(Ap)))
#v_cache = similar(v0)
#L = AffineDiffEqOperator{Float64}(Ap, fp, v_cache)
#function advection(du,u,p,t)
#    du .=  Ap*u .+ fp
#end
L = DiffEqArrayOperator(Ap)

# initialise the plot and define the update function
scene = AbstractPlotting.lines(x, v0)
display(scene)
lineplot = scene[end]
function plotter(integrator)
    v = copy(v0)
    v[iidxs] = integrator.u
    v[bidxs] = gp - Bp*integrator.u
    lineplot[2] = v
    AbstractPlotting.update!(scene)
    sleep(1/60)
end
cb = PeriodicCallback(plotter, Δt)

# define the ODE problem and solve
prob = ODEProblem(L, u0, tspan)
solver = LinearExponential()#Exprb32()#Tsit5()#
u = solve(prob, solver, callback=cb, saveat=tspan)