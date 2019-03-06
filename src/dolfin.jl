# https://github.com/FEniCS/dolfin/blob/df4eea209a85c461e3dcb7c21d4b4c015c46ecdf/dolfin/io/XMLMesh.cpp

using EzXML

mutable struct DolfinMesh
    celltype::String
    gdim::Int
    vertices#::Vector{Tuple{Float64}}
    cells#::Vector{Tuple{UInt}}
end

struct DolfinMeshValueCollection
    typ
    dim
    cell_indices#::Array{Int}
    local_entities#::Array{Int}
    values#::Array{Int}
end

function readdolfinmesh(filename)

    xml_dolfin = readxml(filename)
    xml_mesh = elements(xml_dolfin.root)[1]

    # Get cell type and geometric dimension
    celltype = xml_mesh["celltype"]
    gdim = parse(Int, xml_mesh["dim"])

    # Get topological dimension
    if celltype=="triangle"
        tdim = 3
    elseif celltype=="tetrahedron"
        tdim = 4
    end

    # Get vertices xml node
    xml_vertices = elements(xml_mesh)[1]

    # Get number of vertices and init editor
    num_vertices = parse(UInt, xml_vertices["size"])
    vertices = Vector{NTuple{gdim, Float64}}(undef, num_vertices)

    # Iterate over vertices and add to mesh
    for vertex in elements(xml_vertices)
        i = parse(Int, vertex["index"])
        if gdim==3
            vertices[i+1] = parse.(Float64, (vertex["x"], vertex["y"], vertex["z"]))
        elseif gdim==2
            vertices[i+1] = parse.(Float64, (vertex["x"], vertex["y"]))
        end
    end

    # Get cells node
    xml_cells = elements(xml_mesh)[2]

    # Get number of cells and init editor
    num_cells = parse(UInt, xml_cells["size"])
    cells = Vector{NTuple{tdim, UInt}}(undef, num_cells)

    # Iterate over cells and add to mesh
    for cell in elements(xml_cells)
        i = parse(Int, cell["index"])
        if tdim==4
            cells[i+1] = parse.(UInt, (cell["v0"], cell["v1"], cell["v2"], cell["v3"])) .+ 1
        elseif tdim==3
            cells[i+1] = parse.(UInt, (cell["v0"], cell["v1"], cell["v2"])) .+ 1
        end
    end

    return DolfinMesh(celltype, gdim, vertices, cells)

end

# entities are :vertex, edge, face, facet, cell
# A Vertex is a MeshEntity of topological dimension 0.
# An Edge is a MeshEntity of topological dimension 1.
# A Face is a MeshEntity of topological dimension 2.
# A Facet is a MeshEntity of topological codimension 1.
#   (i.e. the Facet of a Tetrehedron is a Face, of a Triangle is an Edge)
# A Cell is a MeshEntity of topological codimension 0.
#   (i.e. the Cell of a Tetrahedron is a volume, of a Triangle is an area)

function readdolfindomains(filename)

    xml_dolfin = readxml(filename)
    xml_mesh = elements(xml_dolfin.root)[1]

    # Iterate over data
    xml_domains = elements(xml_mesh)[3]
    domains = Vector{DolfinMeshValueCollection}(undef, length(elements(xml_domains)))

    for j = 1:length(elements(xml_domains))

        domain = elements(xml_domains)[j]

        # Check that node is <mesh_value_collection>
        #if domain.name != "mesh_value_collection"
        #    return
        #end

        # Get attributes and init editor
        typ = domain["type"]
        dim = parse(UInt, domain["dim"])
        num_entries = parse(UInt, domain["size"])
        cell_indices = Vector{UInt}(undef, num_entries)
        local_entities = Vector{UInt}(undef, num_entries)
        values = Vector{UInt}(undef, num_entries)

        # iterate over entries
        for i = 1:num_entries
            entry = elements(domain)[i]
            cell_indices[i] = parse(UInt, entry["cell_index"])
            local_entities[i] = parse(UInt, entry["local_entity"])
            values[i] = parse(UInt, entry["value"])
        end
        domains[j] = DolfinMeshValueCollection(typ, dim, cell_indices, local_entities, values)
    end

    return domains

end

# e0 = (v1, v2)
# e1 = (v0, v1)
# e2 = (v0, v1)
# e3 = (v0, v3)
# e4 = (v1, v3)
# e5 = (v2, v3)
#
# f0 = (e0, e4, e5)
# f1 = (e1, e3, e5)
# f2 = (e2, e3, e4)
# f3 = (e0, e1, e2)

# Initialise facet array to be 4x the size of the tet mesh.
# Define facet list: iterate over cells adding unique facets.
# Any time a facet is encountered more than once, flag it as internal.
# Return only the external ones.

#function getboundary(dolfin_mesh::DolfinMesh)
#  // We iterate over all facets in the mesh and check if they are on
#  // the boundary. A facet is on the boundary if it is connected to
#  // exactly one cell.
#  // Any node within a boundary facet is on the boundary.
#end

function getboundary(dolfin_mesh::DolfinMesh)
    num_vertices = length(dolfin_mesh.vertices)
    num_triangles = length(dolfin_mesh.cells)

    # make an array of all edges
    edge_array = Array{Set{Integer}}(undef, 3, num_triangles)
    for i = 1:num_triangles
        tri = dolfin_mesh.cells[i]
        edge1 = Set{UInt}([tri[1], tri[2]])
        edge2 = Set{UInt}([tri[2], tri[3]])
        edge3 = Set{UInt}([tri[3], tri[1]])
        edge_array[1, i] = edge1
        edge_array[2, i] = edge2
        edge_array[3, i] = edge3
    end

    # check whether there are multiple instances
    # really, I just want to compare rows, but a cell shouldn't have copies of an edge anyway, so columns is fine too
    all_edges = edge_array[:]
    unique_edges = unique(all_edges)
    ind_map = indexin(all_edges, unique_edges)
    num_repeats = zeros(length(unique_edges))
    for i = 1:length(all_edges)
        num_repeats[ind_map[i]] += 1
    end
    is_interior = num_repeats .> 1

    # invert and flag boundary points
    is_boundary = falses(num_vertices)
    for i = 1:length(unique_edges)
        if !is_interior[i]
            edge = unique_edges[i]
            is_boundary[collect(edge)] .= true
        end
    end

    return is_boundary

    # cell_edges = num_cells x 3
    # Generate array of all edges
    # Create a unique version
    # ind = indexin(all_edges, unique_edges)
    # num_repeats =

end
