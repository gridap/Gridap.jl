module FESpaces

export FESpace
export ConformingFESpace
export globalnumbering, computelgidvefs
export interpolate
export Assembler
export ConformingAssembler

using Numa: evaluate
using Numa.Helpers
using Numa.RefFEs
using Numa.Polytopes
using Numa.CellValues
using Numa.CellValues.ConstantCellValues
using Numa.CellValues.Wrappers
using Numa.CellValues.Append
using Numa.Geometry
using Numa.FieldValues

using SparseArrays

using Numa.CellValues.Wrappers: CellVectorByComposition
using Numa.CellMaps
using Numa.CellMaps.Operations: CellFieldFromExpand

import Numa: evaluate!

abstract type MeshConformity end
struct ConformingMesh <: MeshConformity end
struct NonConformingMesh <: MeshConformity end

# Abstract types and interfaces

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank) and the number type `E`
"""
abstract type FESpace{D,Z,T,E} end

nf_dofs(::FESpace) = @abstractmethod

cell_eqclass(::FESpace) = @abstractmethod

num_free_dofs(::FESpace) = @abstractmethod

num_fixed_dofs(::FESpace) = @abstractmethod

MeshConformity(::Type{FESpace})::MeshConformity = @abstractmethod
MeshConformity(this::FESpace) = MeshConformity(typeof(this))

function applyconstraints(::FESpace,
	cellvec::CellVector)::Tuple{CellVector,CellVector{Int}}
	@abstractmethod
end

function applyconstraintsrows(::FESpace,
	cellmat::CellMatrix)::Tuple{CellMatrix,CellVector{Int}}
	@abstractmethod
end

function applyconstraintscols(::FESpace,
	cellmat::CellMatrix)::Tuple{CellMatrix,CellVector{Int}}
	@abstractmethod
end

assemblycellgids(::FESpace)::CellVector{Int} = @abstractmethod

_reffes(::FESpace) = @abstractmethod

_dir_tags(::FESpace) = @abstractmethod

_triangulation(::FESpace) = @abstractmethod

_gridgraph(::FESpace) = @abstractmethod

_labels(::FESpace) = @abstractmethod

struct FESpaceWithDirichletData{D,Z,T,E,V<:FESpace{D,Z,T,E}} <: FESpace{D,Z,T,E}
	fesp::V
	dir_data::Vector{Float64}
end

MeshConformity(::Type{FESpaceWithDirichletData{D,Z,T,E,S}})  where {D,Z,T,E,S} = MeshConformity(S)
# MeshConformity(::Type{FESpaceWithDirichletData{D,Z,T,E,S}}) where {D,Z,T,E,S<:NonConformingFESpace{D,Z,T,E}} = NonConformingMesh

for op in (:nf_dofs, :cell_eqclass, :num_free_dofs, :num_fixed_dofs, :dir_tags,
	:_reffes, :_triangulation, :_gridgraph, :_labels)
	@eval begin
		$op(this::FESpaceWithDirichletData) = $op(this.fesp)
	end
end

dir_data(this::FESpaceWithDirichletData) = this.dir_data

function TestFESpace(this::FESpace)
  dv = zeros(Float64,num_fixed_dofs(this))
  return FESpaceWithDirichletData(this, dv)
end

function TrialFESpace( this::FESpace{D}, fun::Vector{Function}, labels::FaceLabels) where {D}
	dv = interpolate_dirichlet_data(fun, this)
	# @santiagobadia : Put labels in FESPace
  return FESpaceWithDirichletData(this, dv)
end

include("ConformingFESpace.jl")
include("FEFunctions.jl")
include("Interpolate.jl")
include("GlobalDofs.jl")
include("Assemblers.jl")

end # module FESpaces
