module FESpaces

export FESpace
export globalnumbering, computelgidvefs

using Numa.Meshes
using Numa.RefFEs
using Numa.Polytopes

# Abstract types and interfaces

# @santiagobadia : To be defined with @fverdugo

# Concrete structs

"""
FE Space structure, where only one RefFE is possible in the whole mesh (to be improved in the future)
"""
struct FESpace
	reffe::LagrangianRefFE
	mesh::Mesh
	l2giddof::Array{Array{Int64,1},1}
	numgdof::Int64
end

# Methods

include("FESpacesMethods.jl")

end # module FESpaces
