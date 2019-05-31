include("MultiFESpaces.jl")

module MultiAssemblers

using Gridap
using Gridap.Helpers
using Gridap.FESpaces
using Gridap.MultiCellArrays
using ..MultiFESpaces
using ..MultiFESpaces: _compute_offsets

using SparseArrays

export MultiAssembler
export MultiSparseMatrixAssembler

import Gridap.Assemblers: assemble
import Gridap.Assemblers: assemble!
import Gridap.Assemblers: SparseMatrixAssembler

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
  testfesps::MultiFESpace{E,S}, trialfesps::MultiFESpace{E,S}) where {E,S}
  @notimplementedif S != ConsequtiveMultiFieldStyle
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

function assemble(this::MultiSparseMatrixAssembler{E},mcm::MultiCellMatrix) where {E}
  V = this.testfesps
  U = this.trialfesps
  mf_vals = apply_constraints_rows(V, mcm)
  mf_rows = celldofids(V)
  mf_vals = apply_constraints_cols(U,mf_vals)
  mf_cols = celldofids(U)
  offsets_row = _compute_offsets(V)
  offsets_col = _compute_offsets(U)
  i_to_fieldid = mf_vals.fieldids
  aux_row = Int[]; aux_col = Int[]; aux_val = E[]
  _assemble_mat!(
    aux_row,aux_col,aux_val,mf_rows,mf_cols,mf_vals,
    i_to_fieldid,offsets_row,offsets_col)
  sparse(aux_row, aux_col, aux_val)
end

function assemble!(
  mat::SparseMatrixCSC{E}, this::MultiSparseMatrixAssembler{E}, vals::MultiCellMatrix) where E
  # This routine can be optimized a lot taking into a count the sparsity graph of mat
  # For the moment we create an intermediate matrix and then transfer the nz values
  m = assemble(this,vals)
  mat.nzval .= m.nzval
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

function _assemble_mat!(
  aux_row,aux_col,aux_val,mf_rows,mf_cols,mf_vals,
  i_to_fieldid,offsets_row,offsets_col)

  for (mf_rows_c,mf_cols_c,mf_vals_c) in zip(mf_rows,mf_cols,mf_vals)
    for (i,vals_c) in enumerate(mf_vals_c)
      ifield,jfield = i_to_fieldid[i]
      rows_c = mf_rows_c[ifield]
      cols_c = mf_cols_c[jfield]
      offset_row = offsets_row[ifield]
      offset_col = offsets_col[jfield]
      for (lidcol,gidcol) in enumerate(cols_c)
        if gidcol > 0
          for (lidrow,gidrow) in enumerate(rows_c)
            if gidrow > 0
              push!(aux_row, gidrow+offset_row)
              push!(aux_col, gidcol+offset_col)
              push!(aux_val, vals_c[lidrow,lidcol])
            end
          end
        end
      end
    end
  end

end

end # module MultiAssemblers
