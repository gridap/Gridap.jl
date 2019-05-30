module Assemblers

using Gridap.FESpaces

using Gridap.Helpers
using Gridap.CellValues

using SparseArrays

export Assembler
export SparseMatrixAssembler
export assemble
export assemble!

#@fverdugo for the moment the abstract interface of Assembler
# (and therefore its concrete implementations)
# assumes single field, and single term

"""
Abstract assembly operator
Parametrized by the type of returned matrix and vector
"""
abstract type Assembler{M<:AbstractMatrix,V<:AbstractVector} end

"""
Assembly of a vector allocating output
"""
function assemble(::Assembler{M,V},::CellVector)::V where {M,V}
  @abstractmethod
end

"""
Assembly of a matrix allocating output
"""
function assemble(::Assembler{M,V},::CellMatrix)::M where {M,V}
  @abstractmethod
end

"""
In-place assembly of a vector (allows optimizations)
"""
function assemble!(::V,::Assembler{M,V},::CellVector)::V where {M,V}
  @abstractmethod
end

"""
In-place assembly of a matrix (allows a LOT of optimizations)
"""
function assemble!(::M,::Assembler{M,V},::CellMatrix)::M where {M,V}
  @abstractmethod
end

"""
Assembler that produces SparseMatrices from the SparseArrays package
"""
struct SparseMatrixAssembler{E} <: Assembler{SparseMatrixCSC{E,Int},Vector{E}}
  testfesp::FESpace
  trialfesp::FESpace
end

function SparseMatrixAssembler(test::FESpace{D,Z,T}, trial::FESpace{D,Z,T}) where {D,Z,T}
  E = eltype(T)
  SparseMatrixAssembler{E}(trial,test)
end

function assemble(this::SparseMatrixAssembler{E}, vals::CellVector) where E
  n = num_free_dofs(this.testfesp)
  vec = zeros(E,n)
  assemble!(vec,this,vals)
  vec
end

function assemble!(
  vec::Vector{E},this::SparseMatrixAssembler{E}, vals::CellVector) where E
  _vals = apply_constraints(this.testfesp, vals)
  _rows = celldofids(this.testfesp)
  vec .= zero(E)
  _assemble_vector!(vec, _vals, _rows)
end

function _assemble_vector!(vec,vals,rows)
  for (rows_c,vals_c) in zip(rows,vals)
    for (i,gid) in enumerate(rows_c)
      if gid > 0
        vec[gid] += vals_c[i]
      end
    end
  end
end

function assemble(this::SparseMatrixAssembler{E}, vals::CellMatrix) where E
  _vals = apply_constraints_rows(this.testfesp, vals)
  rows_m = celldofids(this.testfesp)
  _vals = apply_constraints_cols(this.trialfesp, _vals)
  cols_m = celldofids(this.trialfesp)
  args = _assemble_sparse_matrix_values(_vals,rows_m,cols_m,Int,E)
  sparse(args...)
end

function _assemble_sparse_matrix_values(vals,rows,cols,I,E)
  aux_row = I[]; aux_col = I[]; aux_val = E[]
  for (rows_c, cols_c, vals_c) in zip(rows,cols,vals)
     for (j,gidcol) in enumerate(cols_c)
       if gidcol > 0
        for (i,gidrow) in enumerate(rows_c)
          if gidrow > 0
            push!(aux_row, gidrow)
            push!(aux_col, gidcol)
            push!(aux_val, vals_c[i,j])
          end
        end
      end
    end
  end
  (aux_row, aux_col, aux_val)
end

function assemble!(
  mat::SparseMatrixCSC{E}, this::SparseMatrixAssembler{E}, vals::CellMatrix) where E
  # This routine can be optimized a lot taking into a count the sparsity graph of mat
  # For the moment we create an intermediate matrix and then transfer the nz values
  m = assemble(this,vals)
  mat.nzval .= m.nzval
end

# Draft of multi field assembler
#function _assemble_sparse_matrix_values(mf_vals,mf_rows,mf_cols,I,E)
#  aux_row = I[]; aux_col = I[]; aux_val = E[]
#  for (mf_rows_c, mf_cols_c, mf_vals_c) in zip(mf_rows,mf_cols,mf_vals)
#    for (vals_c, (ifield, jfield)) in eachblock(mf_vals_c)
#      rows_c = mf_rows_c[ifield]
#      cols_c = mf_cols_c[jfield]
#      row_offset = row_offsets[ifield]
#      col_offset = col_offsets[jfield]
#      _asseble_cell_values!(aux_row,aux_col,aux_val,rows_c,cols_c,vals_c,col_offset,row_offset)
#    end
#  end
#  (aux_row, aux_col, aux_val)
#end
#
#function _asseble_cell_values!(aux_row,aux_col,aux_val,rows_c,cols_c,vals_c,col_offset,row_offset)
# for (j,gidcol) in enumerate(cols_c)
#   if gidcol > 0
#     for (i,gidrow) in enumerate(rows_c)
#       if gidrow > 0
#         push!(aux_row, gidrow+row_offset)
#         push!(aux_col, gidcol+col_offset)
#         push!(aux_val, vals_c[i,j])
#       end
#     end
#   end
# end
#end

end # module Assemblers
