module FESpaces

export FESpace
export globalnumbering, computelgidvefs

using Numa.Meshes
using Numa.RefFEs
using Numa.Polytopes

# Abstract types and interfaces

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank) and the number type `E`
"""
abstract type FESpace{D,Z,T,E} end

ndofs(::FESpace) = @abstractmethod

function maprows(::FESpace{D,Z,T,E}, ::CellVector{E},
	::CellValue{Int})::SubAssembledVector{E}  where {D,Z,T,E}
  @abstractmethod
end # instead of maprows, we can call it subassembledvector

# @santiagobadia: How do we express the output? Do we want a new structure that
# will aggregate this information, which is a lazy assembled vector? I would say
# yes. Thus, everything is mutable. Doing assemble(::SubAssembledVector) the
# computation will start

struct SubAssembledVector{E} where E
	cell_values::CellVector{E}
	assembly_op::CellVector{Int}
end

struct SubAssembledMatrix{E} where E
	cell_values::CellMatrix{E}
	assembly_op_rows::CellVector{Int}
	assembly_op_cols::CellVector{Int}
end
# We would need a SubAssembledMatrix constructor that invokes maprows and
# mapcolumns

# @santiagobadia : When defined the vector, these will be straightforward
# function maprows for matrices
# function mapcolums (for matrices only)

# I think, as in many other libraries as Fenics, FEMPAR, etc, we can consider
# the concept of System, that stores both Vector and Matrix. We can define it,
# which will aggregate a SubAssembledVector and a SubAssembledMatrix. In this
# case, we could store assembly_op_rows only once. We would need to define a
# constructor for vector and natrix together.

# So instead of the previous structs, we could have this one

struct SubAssembledSystem{E} where E
	cell_rhs::CellVector{E}
	cell_matrix::CellMatrix{E}
	assembly_op_rows::CellVector{Int}
	assembly_op_cols::CellVector{Int}
end

# A FE Space would be...

struct ConformingFESpace{D,Z,T} <: FESpace{D,Z,T,Float64} where {D,Z,T}
	# For the moment, I am not considering E (to think)
	reffe::LagrangianRefFE{D,T}
	mesh::Mesh{D,Z}
	dof_l_to_g::CellVector{Int}
	numdofs::Int (not sure it can be provided withou)
end

# Old stuff

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
0)
