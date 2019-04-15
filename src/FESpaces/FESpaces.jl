module FESpaces

export FESpace
export globalnumbering, computelgidvefs

using Numa.Helpers
using Numa.RefFEs
using Numa.Polytopes
using Numa.CellValues
using Numa.Geometry

using Numa.Meshes
using SparseArrays

using Numa.CellValues: CellVectorByComposition

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

function ConformingAssembler(fesp::FESpace)
	gldofs = globaldofs(fesp)
	ndofs = gldofs[end][end]
	cell_to_dofs = CellVectorByComposition(fesp.mesh.cellvefs, gldofs)
	ConformingAssembler{Float64}(cell_to_dofs, cell_to_dofs, ndofs)
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

function assemble(this::Assembler, vals::CellVector{T}) where T
  rows_m = this.assembly_op_rows
	# @santiagobadia : Evaluate efficiency, best way to do it in Julia
	# without pre-allocate loop?
  aux_row = []; aux_vals = []
  for vals_c in vals
    aux_vals = [aux_vals..., vals_c...]
  end
  for rows_c in rows_m
    aux_row = [ aux_row..., rows_c...]
  end
  return Array(sparsevec(aux_row, aux_vals))
end

function assemble(this::Assembler, vals::CellMatrix{T}) where T
  rows_m = this.assembly_op_rows
  cols_m = this.assembly_op_cols
	# @santiagobadia : Evaluate efficiency, best way to do it in Julia
	# without pre-allocate loop?
  aux_row = []; aux_col = []; aux_vals = []
  for vals_c in vals
    aux_vals = [aux_vals..., vec(vals_c)...]
  end
  for (rows_c, cols_c) in zip(rows_m,cols_m)
    for I in Iterators.product(rows_c, cols_c)
      aux_row = [aux_row..., I[1]]
      aux_col = [aux_col..., I[2]]
    end
  end
  return sparse(aux_row, aux_col, aux_vals)
end

include("FESpacesMethods.jl")

end # module FESpaces
