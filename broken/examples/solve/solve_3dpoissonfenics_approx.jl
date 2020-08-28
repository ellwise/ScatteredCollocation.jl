# All problems: https://fenicsproject.org/docs/dolfin/1.5.0/python/demo/index.html
# All meshes: https://fenicsproject.org/pub/data/meshes/
# Problem adapted from: https://fenicsproject.org/docs/dolfin/1.5.0/python/demo/documented/bcs/python/documentation.html
# Mesh downloaded from: https://fenicsproject.org/pub/data/meshes/aneurysm.xml

using Plots
using ScatteredCollocation

meshfile = joinpath(@__DIR__,"../../meshes/dolfin_fine.xml")

dolfin_mesh = ScatteredCollocation.readdolfinmesh(meshfile)

#domains = ScatteredCollocation.readdolfindomains(meshfile)

is_boundary = ScatteredCollocation.getboundary(dolfin_mesh)

x_bnd = dolfin_mesh.vertices[is_boundary]
x_int = dolfin_mesh.vertices[.!is_boundary]

plotly()
scatter(map(y->y[1], x_int), map(y->y[2], x_int), markersize=2)#, map(y->y[3], x)
scatter!(map(y->y[1], x_bnd), map(y->y[2], x_bnd), markersize=2)#, map(y->y[3], x)
