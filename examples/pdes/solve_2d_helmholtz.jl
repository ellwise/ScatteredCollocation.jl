using ScatteredCollocation
using Makie
using Arpack
SC = ScatteredCollocation

meshfile = joinpath(@__DIR__,"../../meshes/dolfin_fine.xml")
numeigs = 30
k = 3
accuracy_order = 1
tol = 1e-4
ε = 1e-3
Δt = 0.01
Δx = 0.05
T1 = Float64
T2 = SparseMatrixCSC{T1}

# generate collocation points
#dolfin_mesh = ScatteredCollocation.readdolfinmesh(meshfile)
#x = [[y...] for y in dolfin_mesh.vertices]
#num_points = length(x)
#is_boundary = ScatteredCollocation.getboundary(dolfin_mesh)
x, is_boundary = ScatteredCollocation.griddedSquarePoints(Δx/2); x *= 2
#x, is_boundary = SC.concentricDiscPoints(Δx)
num_points = length(x)
tree = KDTree(hcat(map(y->[y...], x)...))
idxs = collect(1:num_points)
bidxs = idxs[is_boundary]
iidxs = setdiff(idxs, bidxs)

# assemble interior/boundary operators and define RHSs
sca = Assembler(k, accuracy_order, tree)
∇² = derivative(T2, sca, [2, 0], iidxs) + derivative(T2, sca, [0, 2], iidxs)
E  = derivative(T2, sca, [0, 0], bidxs)
f = zeros(T1, length(iidxs))
g = zeros(T1, length(bidxs))

# define initial condition
x0 = [0.75, 0.8]
p0 = exp.(-norm.(x.-Ref(x0),2).^2 ./ (0.05^2)) # pressure

# project the BCs onto the interior
∇² = ∇²[:, vcat(iidxs, bidxs)]
E  =  E[:, vcat(iidxs, bidxs)]
ℒ, ℬ, fp, gp = projbc(∇², E, f, g)
dropzeros!(ℒ)
dropzeros!(ℬ)
fracnz = nnz(ℒ)/prod(size(ℒ))
print("Sparsity = $(1-fracnz)")

p = symrcm(ℒ)
ip = invperm(p)
ℒ = ℒ[p, p]
ℬ = ℬ[:,p]

# compute the eigenvalues and eigensolutions
λ, ϕ = Arpack.eigs(ℒ; nev=numeigs, which=:SM)

# initialise the plot and define the update function
u = Node(zeros(T1, length(idxs)))
clims = maximum(real.(ϕ))*[-1, 1]
scene = AbstractPlotting.Scene(scale_plot=false)
scplot!(scene, tree, u; colorrange=clims, colormap=:pu_or)
display(scene)

for j = 1:numeigs
    unew = zeros(T1, length(idxs))
    unew[iidxs] = real.(ϕ[ip,j])
    unew[bidxs] = gp-ℬ*real.(ϕ[:,j])
    push!(u, unew)
    sleep(0.25)
end