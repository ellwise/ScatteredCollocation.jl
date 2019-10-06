using ScatteredCollocation, Makie

dx = 0.05
xx, bnd = ScatteredCollocation.concentricDiscPoints(dx)
x = collect(range(-1.0,stop=1.0,length=convert(Integer,cld(2,dx))))
#x = [[y[1]-0.5, y[2]-0.5] for y in x]

uu = ScatteredCollocation.phi.(xx, 3, Ref([0, 0]), false, Ref([1, 0]))

scene = AbstractPlotting.Scene()
scatter!(scene, map(y->y[1], xx), map(y->y[2], xx), uu, markersize = 0.03)
display(scene)