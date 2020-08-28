using ScatteredCollocation, Makie
const SC = ScatteredCollocation

Δx = 0.05
x, bnd = SC.concentricDiscPoints(Δx)

maxorder = 3
numdims = 2
npoly = SC.countmonomials(maxorder; numdims=numdims)
m = SC.MonomialExponents{numdims,Int}(npoly)
poly_exponents = collect(Iterators.take(m, npoly))

islog = false
diff_orders = [0, 2]
ks = (0:4)# .+ sum(diff_orders)

scenes = []
for k in ks
    row = []
    for ps in poly_exponents
        u = SC.phi.(x, k, Ref(ps), islog, Ref(diff_orders))
        clims = 1.0*maximum(abs.(u))*[-1, 1]
        scene = AbstractPlotting.Scene(scale_plot=false)
        scatter!(scene, map(y->y[1], x), map(y->y[2], x), markersize = 1.1*Δx, color=u,
            colorrange=clims, colormap=:pu_or, show_axis=false)
        push!(row, scene)
    end
    push!(scenes, row)
end
tot = hbox([vbox([scene for scene in row]) for row in scenes])
display(tot)
