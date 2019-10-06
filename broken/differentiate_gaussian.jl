using ScatteredCollocation
using Makie
using LinearAlgebra

meshfile = joinpath(@__DIR__,"../../meshes/dolfin_fine.xml")
k = 3
accuracy_order = 3
Δt = 0.02
dx = 0.05
T1 = Float64
T2 = SparseMatrixCSC{T1}

# generate collocation points
#dolfin_mesh = ScatteredCollocation.readdolfinmesh(meshfile)
#x = [[y...] for y in dolfin_mesh.vertices]
#num_points = length(x)
#is_boundary = ScatteredCollocation.getboundary(dolfin_mesh)
#x, is_boundary = ScatteredCollocation.griddedSquarePoints(0.05)
x, is_boundary = ScatteredCollocation.concentricDiscPoints(dx)
num_points = length(x)
tree = KDTree(hcat(map(y->[y...], x)...),reorder=false)
idxs = collect(1:num_points)
bidxs = idxs[is_boundary]
iidxs = setdiff(idxs, bidxs)

# assemble operators
sca = SCAssembler3(k, accuracy_order, tree)
Dx = derivative(T2, sca, [1, 0], idxs)
#scene = AbstractPlotting.Scene()
#heatmap!(scene, Matrix(abs.(Dx)))
#display(scene)
#dfaasdfsfda
Dy = derivative(T2, sca, [0, 1], idxs)
E  = derivative(T2, sca, [0, 0], idxs)

# define initial condition
x0 = [0.3, 0.3]
p0 = exp.(-norm.(x.-Ref(x0),2).^2 ./ (0.1^2)) # pressure

markersize = 0.01*Dy*p0
scene = AbstractPlotting.Scene()
scatter!(map(y->y[1], x), map(y->y[2], x), markersize=markersize, color=:blue)
display(scene)




afsafdfads
# project the BCs onto the interior
A = [A[:,idxs] A[:,idxs.+num_points] A[:,iidxs.+2*num_points] A[:,bidxs.+2*num_points]]
B = [B[:,idxs] B[:,idxs.+num_points] B[:,iidxs.+2*num_points] B[:,bidxs.+2*num_points]]
Ap, Bp, fp, gp = projbc(A, B, f, g)
dropzeros!(Ap)
dropzeros!(Bp)
print("Sparsity = "); println(1-nnz(Ap)/prod(size(Ap)))
L = DiffEqArrayOperator(Ap)

# initialise the plot and define the update function
pint = Node(p0[iidxs])
pbnd = Node(p0[bidxs])
markersizeint = lift(x->abs.(0.1*x), pint)
markersizebnd = lift(x->abs.(0.1*x), pbnd)
scene = AbstractPlotting.Scene()
scatter!(map(y->y[1], x[iidxs]), map(y->y[2], x[iidxs]), markersize=markersizeint, color=:blue)
scatter!(map(y->y[1], x[bidxs]), map(y->y[2], x[bidxs]), markersize=markersizebnd, color=:red)
display(scene)
function plotter(integrator)
    println(Δt)
    push!(pint, integrator.u[2*num_points+1:end])
    push!(pbnd, gp - Bp*integrator.u)
    sleep(1/60)
end
cb = PeriodicCallback(plotter, Δt)

# define the ODE problem and solve
prob = ODEProblem(L, u0, tspan)
solver = KenCarp4()
u = solve(prob, solver, callback=cb, saveat=tspan)