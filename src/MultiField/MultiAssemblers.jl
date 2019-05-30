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

function assemble(this::MultiSparseMatrixAssembler{E},vals::MultiCellVector) where E
  n = num_free_dofs(this.testfesps)
  vec = zeros(E,n)
  assemble!(vec,this,vals)
  vec
end

function assemble!(vec::Vector{E},this::MultiSparseMatrixAssembler{E},mcv::MultiCellVector) where E
  vec .= zero(E)
  V = this.testfesps
  mf_vals = apply_constraints(V,mcv)
  mf_rows = celldofids(V)
  offsets = _compute_offsets(V)
  i_to_fieldid = mf_vals.fieldids
  _assemble_vec!(vec,mf_rows,mf_vals,i_to_fieldid,offsets)
end

function _assemble_vec!(vec,mf_rows,mf_vals,i_to_fieldid,offsets)
  for (mf_rows_c,mf_vals_c) in zip(mf_rows,mf_vals)
    for (i,vals_c) in enumerate(mf_vals_c)
      ifield, = i_to_fieldid[i]
      rows_c = mf_rows_c[ifield]
      offset = offsets[ifield]
      for (lid,gid) in enumerate(rows_c)
        if gid > 0
          vec[gid+offset] += vals_c[lid]
        end
      end
    end
  end
end

function _compute_offsets(U)
  n = length(U)
  offsets = zeros(Int,n)
  for i in 1:(n-1)
    Ui = U[i]
    offsets[i+1] = offsets[i] + num_free_dofs(Ui)
  end
  offsets
end

end # module MultiAssemblers
