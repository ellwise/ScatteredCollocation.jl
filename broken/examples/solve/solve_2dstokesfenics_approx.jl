# All problems: https://fenicsproject.org/docs/dolfin/1.5.0/python/demo/index.html
# All meshes: https://fenicsproject.org/pub/data/meshes/
# Problem adapted from: https://fenicsproject.org/docs/dolfin/1.5.0/python/demo/documented/stokes-taylor-hood/python/documentation.html
# Mesh downloaded from: https://fenicsproject.org/pub/data/meshes/dolfin_fine.xml

using EzXML
using Plots

xdolfin = readxml(joinpath(@__DIR__,"../../meshes/dolfin_fine.xml"))


xmesh = elements(xdolfin.root)[1]
xvertices, xcells = elements(xmesh)

vertices = elements(xvertices)
x = Array{Tuple{Float64, Float64}}(undef, length(vertices))
for vertex in vertices
    i = parse(Int, vertex["index"])
    x[i+1] = parse.(Float64, (vertex["x"], vertex["y"]))
end

plotly()
scatter(map(y->y[1], x), map(y->y[2], x), markersize=2)
