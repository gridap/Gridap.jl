module CLagrangianFESpaces

using Gridap
using Gridap.DOFBases: _length
using Gridap.ConformingFESpaces: _CellField

struct CLagrangianFESpace{D,Z,T} <: FESpace{D,Z,T}
  fnode_to_coord::Vector{Point{D,Float64}}
  dnode_to_coord::Vector{Point{D,Float64}}
end

num_free_dofs(fesp::FESpace{D,Z,L}) = length(fesp.fnode_to_coord)

num_diri_dofs(fesp::FESpace) = @abstractmethod

diri_tags(fesp::FESpace) = @abstractmethod

function apply_constraints(::FESpace, cellvec::CellVector, cellids::CellNumber)::CellVector
  @abstractmethod
end

function apply_constraints_rows(::FESpace, cellmat::CellMatrix, cellids::CellNumber)::CellMatrix
  @abstractmethod
end

function apply_constraints_cols(::FESpace, cellmat::CellMatrix, cellids::CellNumber)::CellMatrix
  @abstractmethod
end

function celldofids(::FESpace)::CellVector{Int}
  @abstractmethod
end

function interpolate_diri_values(::FESpace, funs::Vector{<:Function})::Vector{E} where E
  @abstractmethod
end

function CellField(
  ::FESpace{D,Z,T},free_dofs::AbstractVector{E},diri_dofs::AbstractVector{E})::CellField{Z,T} where {D,Z,T,E}
  @abstractmethod
end

function CellBasis(::FESpace{D,Z,T})::CellBasis{Z,T} where {D,Z,T}
  @abstractmethod
end

end # module
