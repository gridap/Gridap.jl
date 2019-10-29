module MultiAssemblers

using Gridap
using Gridap.Helpers
using Gridap.MultiFESpaces: _compute_offsets

using SparseArrays
using SparseMatricesCSR

export MultiAssembler
export MultiSparseMatrixAssembler

import Gridap.Assemblers: assemble
import Gridap.Assemblers: assemble!
import Gridap.Assemblers: SparseMatrixAssembler
import Gridap.Assemblers: SparseMatrixAssembler
using  Gridap.Assemblers: sparse_from_coo

abstract type MultiAssembler{M<:AbstractMatrix,V<:AbstractVector} end

function assemble(
  ::MultiAssembler{M,V},
  ::Vararg{Tuple{<:MultiCellVector,<:CellNumber}})::V where {M,V}
  @abstractmethod
end

function assemble(
  ::MultiAssembler{M,V},
  ::Vararg{Tuple{<:MultiCellMatrix,<:CellNumber,<:CellNumber}})::M where {M,V}
  @abstractmethod
end

function assemble!(
  ::V,
  ::MultiAssembler{M,V},
  ::Vararg{Tuple{<:MultiCellVector,<:CellNumber}}) where {M,V}
  @abstractmethod
end

function assemble!(
  ::M,
  ::MultiAssembler{M,V},
  ::Vararg{Tuple{<:MultiCellMatrix,<:CellNumber,<:CellNumber}}) where {M,V}
  @abstractmethod
end

function assemble(
  a::MultiAssembler,cv::MultiCellVector)
  l = length(cv)
  ide = IdentityCellNumber(Int,l)
  assemble(a,(cv,ide))
end

function assemble(
  a::MultiAssembler,cv::MultiCellMatrix)
  l = length(cv)
  ide = IdentityCellNumber(Int,l)
  assemble(a,(cv,ide,ide))
end

function assemble!(r,a::MultiAssembler,cv::MultiCellVector)
  l = length(cv)
  ide = IdentityCellNumber(Int,l)
  assemble!(r,a,(cv,ide))
end

function assemble!(r,a::MultiAssembler,cv::MultiCellMatrix)
  l = length(cv)
  ide = IdentityCellNumber(Int,l)
  assemble!(r,a,(cv,ide,ide))
end

"""
MultiAssembler that produces SparseMatrices from the SparseArrays package
"""
struct MultiSparseMatrixAssembler{E,M} <: MultiAssembler{M,Vector{E}}
  testfesps::MultiFESpace{E}
  trialfesps::MultiFESpace{E}

  function MultiSparseMatrixAssembler(
    ::Type{M},
    testfesp::MultiFESpace{E},
    trialfesp::MultiFESpace{E}) where {M<:AbstractSparseMatrix,E}
    new{E,M}(testfesp,trialfesp) 
  end
end

function MultiSparseMatrixAssembler(
  testfesp::MultiFESpace{E},
  trialfesp::MultiFESpace{E}) where {E}
  MultiSparseMatrixAssembler(SparseMatrixCSC, testfesps, trialfesp)
end

function MultiSparseMatrixAssembler(
  ::Type{M},
  testfesps::Vector{<:FESpaceWithDirichletData},
  trialfesps::Vector{<:FESpaceWithDirichletData}) where {M<:AbstractSparseMatrix}
  V = MultiFESpace(testfesps)
  U = MultiFESpace(trialfesps)
  MultiSparseMatrixAssembler(M,V,U)
end

function MultiSparseMatrixAssembler(
  testfesps::Vector{<:FESpaceWithDirichletData},
  trialfesps::Vector{<:FESpaceWithDirichletData})
  MultiSparseMatrixAssembler(SparseMatrixCSC,V,U)
end

function SparseMatrixAssembler(
  ::Type{M},
  testfesps::Vector{<:FESpaceWithDirichletData},
  trialfesps::Vector{<:FESpaceWithDirichletData}) where {M<:AbstractSparseMatrix}
  MultiSparseMatrixAssembler(M,testfesps,trialfesps)
end

function SparseMatrixAssembler(
  testfesps::Vector{<:FESpaceWithDirichletData},
  trialfesps::Vector{<:FESpaceWithDirichletData})
  MultiSparseMatrixAssembler(SparseMatrixCSC,testfesps,trialfesps)
end

function MultiSparseMatrixAssembler(
  ::Type{M},
  testfesps::MultiFESpace{E,S}, 
  trialfesps::MultiFESpace{E,S}) where {M<:AbstractSparseMatrix,E,S}
  @notimplementedif S != ConsequtiveMultiFieldStyle
  MultiSparseMatrixAssembler{E,M}(M,testfesps,trialfesps)
end

function MultiSparseMatrixAssembler(
  testfesps::MultiFESpace{E,S}, 
  trialfesps::MultiFESpace{E,S}) where {E,S}
  @notimplementedif S != ConsequtiveMultiFieldStyle
  MultiSparseMatrixAssembler(SparseMatrixCSC,testfesps,trialfesps)
end

function assemble(
  this::MultiSparseMatrixAssembler{E},
  allvals::Vararg{Tuple{<:MultiCellVector,<:CellNumber}}) where {E}

  n = num_free_dofs(this.testfesps)
  vec = zeros(E,n)
  assemble!(vec,this,allvals...)
  vec
end

function assemble!(
  vec::Vector{E},
  this::MultiSparseMatrixAssembler{E},
  allmcv::Vararg{Tuple{<:MultiCellVector,<:CellNumber}}) where {E}

  vec .= zero(E)
  V = this.testfesps
  offsets = _compute_offsets(V)
  _mf_rows = celldofids(V)
  for (mcv,cellids) in allmcv
    mf_vals = apply_constraints(V,mcv,cellids)
    i_to_fieldid = mf_vals.fieldids
    mf_rows = reindex(_mf_rows,cellids)
    _assemble_vec!(vec,mf_rows,mf_vals,i_to_fieldid,offsets)
  end
end

function assemble(
  this::MultiSparseMatrixAssembler{E,M},
  allmcm::Vararg{Tuple{<:MultiCellMatrix,<:CellNumber,<:CellNumber}}) where {E,M}

  V = this.testfesps
  U = this.trialfesps
  offsets_row = _compute_offsets(V)
  offsets_col = _compute_offsets(U)
  _mf_rows = celldofids(V)
  _mf_cols = celldofids(U)
  aux_row = Int[]; aux_col = Int[]; aux_val = E[]

  for (mcm,cellids_row,cellids_col) in allmcm
    mf_vals = apply_constraints_rows(V, mcm, cellids_row)
    mf_rows = reindex(_mf_rows, cellids_row)
    mf_vals = apply_constraints_cols(U,mf_vals,cellids_col)
    mf_cols = reindex(_mf_cols, cellids_col)
    i_to_fieldid = mf_vals.fieldids
    _assemble_mat!(M,
      aux_row,aux_col,aux_val,mf_rows,mf_cols,mf_vals,
      i_to_fieldid,offsets_row,offsets_col)
  end
  num_rows = num_free_dofs(this.testfesps)
  num_cols = num_free_dofs(this.trialfesps)
  finalize_coo!(M,aux_row,aux_col,aux_val,num_rows,num_cols)
  sparse_from_coo(M,aux_row, aux_col, aux_val)
end

function assemble!(
  mat::AbstractSparseMatrix{E},
  this::MultiSparseMatrixAssembler{E},
  allvals::Vararg{Tuple{<:MultiCellMatrix,<:CellNumber,<:CellNumber}}) where E
  # This routine can be optimized a lot taking into a count the sparsity graph of mat
  # For the moment we create an intermediate matrix and then transfer the nz values
  m = assemble(this,allvals...)
  nonzeros(mat) .= nonzeros(m)
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

function _assemble_mat!(::Type{M},
  aux_row,aux_col,aux_val,mf_rows,mf_cols,mf_vals,
  i_to_fieldid,offsets_row,offsets_col) where {M<:AbstractSparseMatrix}

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
              push_coo!(M,aux_row, aux_col, aux_val, 
                    gidrow+offset_row, gidcol+offset_col, vals_c[lidrow,lidcol])
            end
          end
        end
      end
    end
  end

end

end # module MultiAssemblers
