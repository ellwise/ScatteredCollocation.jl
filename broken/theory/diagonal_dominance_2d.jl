using ScatteredCollocation
using Makie
using LinearAlgebra
using SparseArrays, NearestNeighbors
SC = ScatteredCollocation

k = 3
accuracy_order = 3
buffer = 0
Δx = 0.1/2
T1 = Float64
T2 = SparseMatrixCSC{T1}

x, is_boundary = SC.concentricDiscPoints(Δx)
x, is_boundary = SC.griddedSquarePoints(Δx)
num_points = length(x)
tree = KDTree(hcat(map(y->[y...], x)...))
idxs = 1:num_points

# assemble interior/boundary operators and define RHSs
sca = Assembler(k, accuracy_order, tree)
∂ₓ = derivative(T2, sca, [0, 1], idxs)
∇² = derivative(T2, sca, [2, 0], idxs) + derivative(T2, sca, [0, 2], idxs)

p0 = Node(zeros(num_points))
clims = 0.1*norm(∂ₓ[:], 2)*[-1, 1]
scene = AbstractPlotting.Scene(scale_plot=false)
scatter!(map(y->y[1], x), map(y->y[2], x), markersize=Δx, color=p0, colorrange=clims, colormap=:pu_or)
display(scene)

for i in 1:num_points

    u = Vector(∂ₓ[i, :])
    #u[iszero.(u)] .= NaN
    push!(p0, u)
    
    sleep(1/60)
end

