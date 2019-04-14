module FESpaces

export FESpace
export globalnumbering, computelgidvefs

using Numa.Helpers
using Numa.RefFEs
using Numa.Polytopes
using Numa.CellValues
using Numa.Geometry

using Numa.Meshes

# Abstract types and interfaces

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank) and the number type `E`
"""
abstract type FESpace{D,Z,T,E} end

function globaldofs(::FESpace{D,Z,T,E})::CellVector{Int}  where {D,Z,T,E}
  @abstractmethod
end

"""
Abstract assembly operator
"""
abstract type Assembler{E} end

ndofs(::FESpace) = @abstractmethod

function assemblevector(::CellVector{E})::Vector{E} where {E}
	@abstractmethod
end

function assemblematrix(::CellMatrix{E})::Matrix{E} where {E}
	@abstractmethod
end

struct ConformingAssembler{E} <: Assembler{E}
	assembly_op_rows::CellVector{Int}
	assembly_op_cols::CellVector{Int}
	num_dofs::Int
end

"""
Conforming FE Space, where only one RefFE is possible in the whole mesh
"""
struct ConformingFESpace{D,Z,T} <: FESpace{D,Z,T,Float64}
	# For the moment, I am not considering E (to think)
	reffe::LagrangianRefFE{D,T}
	mesh::Mesh{D}
end

# Methods



include("FESpacesMethods.jl")

end # module FESpaces
