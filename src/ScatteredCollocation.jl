#__precompile__()

"""
Scattered collocation provides methods for generating inner product weights
corresponding to differential operators on scattered stencil points.
"""
module ScatteredCollocation

using AMD: amd

using SparseArrays
#using LinearAlgebra
using LinearAlgebra: pinv, norm, Diagonal, Symmetric, I
using IterativeSolvers: gmres
#using NearestNeighbors
using NearestNeighbors: NNTree, knn
#using StaticArrays
using StaticArrays: MVector, SVector
#using Makie
using Makie: lines!, scatter!, heatmap!
using IncompleteLU: ilu
#using AlgebraicMultigrid: ruge_stuben, aspreconditioner
#using Preconditioners
##using Arpack: svds

include("./geometry.jl")
include("./basis.jl")
include("./utils.jl")
include("./assembly.jl")
#include("./dolfin.jl")
include("./CuthillMcKee.jl")
include("./experimental.jl")
include("./sparsealgebra.jl")

export AffineOperator
export derivative, projbc, scplot!, scspy!, symrcm, fnz, addop!, build, projbc2

end # module
