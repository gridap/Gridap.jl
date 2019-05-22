module Assemblers

using Gridap.FESpaces

using Gridap.Helpers
using Gridap.CellValues

using SparseArrays

export Assembler
export SparseMatrixAssembler
export assemble

#@fverdugo for the moment the abstract interface of Assembler
# (and therefore its concrete implementations)
# assumes single field, and single term

"""
Abstract assembly operator
Parametrized by the type of returned matrix and vector
"""
abstract type Assembler{M<:AbstractMatrix,V<:AbstractVector} end

function assemble(::Assembler{M,V},::CellVector)::V where {M,V}
  @abstractmethod
end

function assemble(::Assembler{M,V},::CellMatrix)::M where {M,V}
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

function assemble(this::SparseMatrixAssembler, vals::CellVector)
  _vals, rows_m = apply_constraints(this.testfesp, vals)
  # @santiagobadia : Evaluate efficiency, best way to do it in Julia
  # without pre-allocate loop?
  aux_row = Int[]; aux_vals = Int[]
  for (rows_c,vals_c) in zip(rows_m,_vals)
    for (i,gid) in enumerate(rows_c)
      if gid > 0
        aux_vals = [aux_vals..., vals_c[i]]
        aux_row = [ aux_row..., rows_c[i]...]
      end
    end
  end
  return Array(sparsevec(aux_row, aux_vals))
end

function assemble(this::SparseMatrixAssembler, vals::CellMatrix)
  _vals, rows_m = apply_constraints_rows(this.testfesp, vals)
  _vals, cols_m = apply_constraints_cols(this.trialfesp, _vals, )
  # @santiagobadia : Evaluate efficiency, best way to do it in Julia
  # without pre-allocate loop?
  aux_row = Int[]; aux_col = Int[]; aux_vals = Int[]
  for vals_c in _vals
  end
  for (rows_c, cols_c, vals_c) in zip(rows_m,cols_m,_vals)
    for (i,gidrow) in enumerate(rows_c)
      if gidrow > 0
        for (j,gidcol) in enumerate(cols_c)
          if gidcol > 0
            aux_row = [aux_row..., rows_c[i]]
            aux_col = [aux_col..., cols_c[j]]
            aux_vals = [aux_vals..., vals_c[i,j]]
          end
        end
      end
    end
  end
  return sparse(aux_row, aux_col, aux_vals)
end

end # module Assemblers
