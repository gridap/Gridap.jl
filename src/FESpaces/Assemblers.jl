module Assemblers

using Gridap

using Gridap.Helpers

using SparseArrays

export Assembler
export SparseMatrixAssembler
export assemble
export assemble!

"""
Abstract assembly operator
Parametrized by the type of returned matrix and vector
"""
abstract type Assembler{M<:AbstractMatrix,V<:AbstractVector} end

"""
Assembly of a vector allocating output
"""
function assemble(
  ::Assembler{M,V},
  ::Vararg{Tuple{<:CellVector,<:CellNumber}})::V where {M,V}
  @abstractmethod
end

"""
Assembly of a matrix allocating output
"""
function assemble(
  ::Assembler{M,V},
  ::Vararg{Tuple{<:CellMatrix,<:CellNumber}})::M where {M,V}
  @abstractmethod
end

"""
In-place assembly of a vector (allows optimizations)
"""
function assemble!(
  ::V,
  ::Assembler{M,V},
  ::Vararg{Tuple{<:CellVector,<:CellNumber}})::V where {M,V}
  @abstractmethod
end

"""
In-place assembly of a matrix (allows a LOT of optimizations)
"""
function assemble!(
  ::M,
  ::Assembler{M,V},
  ::Vararg{Tuple{<:CellMatrix,<:CellNumber}})::M where {M,V}
  @abstractmethod
end

function assemble(a::Assembler,cv::CellArray)
  l = length(cv)
  ide = IdentityCellNumber(Int,l)
  assemble(a,(cv,ide))
end

function assemble!(r,a::Assembler,cv::CellArray)
  l = length(cv)
  ide = IdentityCellNumber(Int,l)
  assemble!(r,a,(cv,ide))
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

function assemble(
  this::SparseMatrixAssembler{E},
  vals::Vararg{Tuple{<:CellVector,<:CellNumber}}) where E

  n = num_free_dofs(this.testfesp)
  vec = zeros(E,n)
  assemble!(vec,this,vals...)
  vec
end

function assemble!(
  vec::Vector{E},
  this::SparseMatrixAssembler{E},
  allvals::Vararg{Tuple{<:CellVector,<:CellNumber}}) where E

  vec .= zero(E)
  rows = celldofids(this.testfesp)
  for (vals, cellids) in allvals
    _vals = apply_constraints(this.testfesp, vals, cellids)
    _rows = reindex(rows,cellids)
    _assemble_vector!(vec, _vals, _rows)
  end
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

function assemble(
  this::SparseMatrixAssembler{E},
  allvals::Vararg{Tuple{<:CellMatrix,<:CellNumber}}) where E

  I = Int
  aux_row = I[]; aux_col = I[]; aux_val = E[]

  _rows_m = celldofids(this.testfesp)
  _cols_m = celldofids(this.trialfesp)

  for (vals,cellids) in allvals
    _vals = apply_constraints_rows(this.testfesp, vals, cellids)
    rows_m = reindex(_rows_m, cellids)
    vals_m = apply_constraints_cols(this.trialfesp, _vals, cellids)
    cols_m = reindex(_cols_m, cellids)
    _assemble_sparse_matrix_values!(
      aux_row,aux_col,aux_val,vals_m,rows_m,cols_m)
  end
  sparse(aux_row,aux_col,aux_val)
end

function _assemble_sparse_matrix_values!(aux_row,aux_col,aux_val,vals,rows,cols)
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
end

function assemble!(
  mat::SparseMatrixCSC{E},
  this::SparseMatrixAssembler{E},
  vals::Vararg{Tuple{<:CellMatrix,<:CellNumber}}) where E
  # This routine can be optimized a lot taking into a count the sparsity graph of mat
  # For the moment we create an intermediate matrix and then transfer the nz values
  m = assemble(this,vals...)
  mat.nzval .= m.nzval
end

end # module
