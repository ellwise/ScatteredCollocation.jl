using ScatteredCollocation, NearestNeighbors
SC = ScatteredCollocation
using LinearAlgebra, Makie
using SparseArrays

numdims = 1
num_nodes = 100
k = 3
T = SparseMatrixCSC{Float64}#Array{Float64}
num_dims = 1
accorder = 2
difforder = 2
maxorder = accorder+difforder
num_poly = SC.countmonomials(maxorder; numdims=numdims)
stencil_size = 2*num_poly
buffer = 0

x = Float64.(collect(1:num_nodes))# .+ 0.1*rand(1)
tree = KDTree(hcat(map(y->[y...], x)...))
sca = Assembler(k, accorder, tree)

idxs = collect(1:num_nodes)
bidxs = [1, num_nodes]
iidxs = setdiff(idxs, bidxs)

A = derivative(T, sca, [difforder], iidxs)
B = derivative(T, sca, [0], bidxs)
A = A[:, vcat(iidxs, bidxs)]
B = B[:, vcat(iidxs, bidxs)]
f = zeros(length(iidxs))
g = zeros(length(bidxs))
Ap, Bp, fp, gp = projbc(A, B, f, g)

#e = eigvals(Dx-gamma*inv(A))
e = eigvals(Array(Ap))
println(maximum(real.(e)))
scene = AbstractPlotting.Scene()#scale_plot=false)
scatter!(real.(e), imag.(e))#, xlim=(-0.015,0.015), ylim=(-2,2))
display(scene)