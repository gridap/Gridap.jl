include("MultiFESpaces.jl")

module MultiAssemblers

using Gridap
using Gridap.Helpers
using Gridap.FESpaces
using Gridap.MultiCellArrays
using ..MultiFESpaces

using SparseArrays

export MultiAssembler
export restrict_to_field
export MultiSparseMatrixAssembler

import Gridap.Assemblers: assemble
import Gridap.Assemblers: assemble!
import Gridap.Assemblers: SparseMatrixAssembler
import Base: length

abstract type MultiAssembler{M<:AbstractMatrix,V<:AbstractVector} end

function assemble(::MultiAssembler{M,V},::MultiCellVector)::V where {M,V}
  @abstractmethod
end

function assemble(::MultiAssembler{M,V},::MultiCellMatrix)::M where {M,V}
  @abstractmethod
end

function assemble!(::V,::MultiAssembler{M,V},::MultiCellVector)::V where {M,V}
  @abstractmethod
end

function assemble!(::M,::MultiAssembler{M,V},::MultiCellMatrix)::M where {M,V}
  @abstractmethod
end

length(::MultiAssembler) = @abstractmethod

function restrict_to_field(::MultiAssembler{M,V},::V,field::Integer)::AbstractVector where {M,V}
  @abstractmethod
end

"""
Assembler that produces SparseMatrices from the SparseArrays package
"""
struct MultiSparseMatrixAssembler{E} <: MultiAssembler{SparseMatrixCSC{E,Int},Vector{E}}
  testfesps::MultiFESpace{E}
  trialfesps::MultiFESpace{E}
end

function MultiSparseMatrixAssembler(
  testfesps::Vector{<:FESpaceWithDirichletData},
  trialfesps::Vector{<:FESpaceWithDirichletData})
  V = MultiFESpace(testfesps)
  U = MultiFESpace(trialfesps)
  MultiSparseMatrixAssembler(V,U)
end

function SparseMatrixAssembler(
  testfesps::Vector{<:FESpaceWithDirichletData},
  trialfesps::Vector{<:FESpaceWithDirichletData})
  MultiSparseMatrixAssembler(testfesps,trialfesps)
end

function MultiSparseMatrixAssembler(
  testfesps::MultiFESpace{E}, trialfesps::MultiFESpace{E}) where E
  @assert length(testfesps) == length(trialfesps)
  MultiSparseMatrixAssembler{E}(testfesps,trialfesps)
end

function assemble(this::MultiSparseMatrixAssembler{E},::MultiCellVector) where E
  n = num_free_dofs(this.testfesps)
  vec = zeros(E,n)
  assemble!(vec,this,vals)
  vec
end

#function assemble!(::Vector{E},::MultiSparseMatrixAssembler{E},::MultiCellVector) where E
#  vec .= zero(E)
#
#
#end

end # module MultiAssemblers
