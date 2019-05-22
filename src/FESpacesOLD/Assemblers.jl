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
  testfesp::FESpace
  trialfesp::FESpace
  num_dofs::Int
end

function ConformingAssembler(this::FESpace)
  ConformingAssembler{Int}(this, this, num_free_dofs(this))
end

function ConformingAssembler(test::FESpace, trial::FESpace)
  @assert trial.num_free_dofs == test.num_free_dofs
  ConformingAssembler{Int}(trial, test, num_free_dofs(trial))
end

# Methods

function assemble(this::Assembler, vals::CellVector{T}) where T
  _vals, rows_m = applyconstraints(this.testfesp, vals)
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

function assemble(this::Assembler, vals::CellMatrix{T}) where T
  _vals, rows_m = applyconstraintsrows(this.testfesp, vals)
  _vals, cols_m = applyconstraintscols(this.trialfesp, _vals, )
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
