
function collect_cell_matrix(a::CollectionOfCellContribution)
  w = []
  r = []
  for cell_contribution in values(a.dict)
    push!(w,get_array(cell_contribution))
    push!(r,get_cell_id(cell_contribution))
  end
  (w,r,r)
end

function collect_cell_matrix(a::CellContribution)
  w = [get_array(a)]
  r = [get_cell_id(a)]
  (w,r,r)
end

function collect_cell_vector(a::CollectionOfCellContribution)
  w = []
  r = []
  for cell_contribution in values(a.dict)
    push!(w,get_array(cell_contribution))
    push!(r,get_cell_id(cell_contribution))
  end
  (w,r)
end

function collect_cell_vector(a::CellContribution)
  w = [get_array(a)]
  r = [get_cell_id(a)]
  (w,r)
end

function collect_cell_matrix_and_vector(
  biform::Union{CellContribution,CollectionOfCellContribution},
  liform::Union{CellContribution,CollectionOfCellContribution})

  collect_cell_matrix_and_vector(
    CollectionOfCellContribution(biform),
    CollectionOfCellContribution(liform))
end

function collect_cell_matrix_and_vector(
  biform::Union{CellContribution,CollectionOfCellContribution},
  liform::Union{CellContribution,CollectionOfCellContribution},
  uhd)

  collect_cell_matrix_and_vector(
    CollectionOfCellContribution(biform),
    CollectionOfCellContribution(liform),
    uhd)
end

function collect_cell_matrix_and_vector(
  biform::CollectionOfCellContribution,liform::CollectionOfCellContribution)

  matvec, mat, vec = _pair_contribution_when_possible(biform,liform)
  matvecdata = collect_cell_matrix(matvec)
  matdata = collect_cell_matrix(mat)
  vecdata = collect_cell_vector(vec)
  (matvecdata, matdata, vecdata)
end

function collect_cell_matrix_and_vector(
  biform::CollectionOfCellContribution,liform::CollectionOfCellContribution,uhd)

  @assert is_a_fe_function(uhd)
  matvec, mat, vec = _pair_contribution_when_possible(biform,liform,uhd)
  matvecdata = collect_cell_matrix(matvec)
  matdata = collect_cell_matrix(mat)
  vecdata = collect_cell_vector(vec)
  (matvecdata, matdata, vecdata)
end

function _pair_contribution_when_possible(biform,liform)
  matvec = Dict{UInt,CellContribution}()
  mat = Dict{UInt,CellContribution}()
  vec = Dict{UInt,CellContribution}()
  for (id,t) in biform.dict
    if haskey(liform.dict,id)
      a = pair_arrays(get_array(biform.dict[id]),get_array(liform.dict[id]))
      matvec[id] = CellContribution(a,t.cell_id,t.object)
    else
      mat[id] = t
    end
  end
  for (id,t) in liform.dict
    if ! haskey(biform.dict,id)
      vec[id] = t
    end
  end
  map(CollectionOfCellContribution, (matvec, mat, vec))
end

function _pair_contribution_when_possible(biform,liform,uhd)
  matvec = Dict{UInt,CellContribution}()
  mat = Dict{UInt,CellContribution}()
  _matvec, _mat, _vec = _pair_contribution_when_possible(biform,liform)
  for (id,t) in _matvec.dict
    cellvals = get_cell_values(uhd,t.cell_id)
    a = attach_dirichlet(get_array(t),cellvals)
    matvec[id] = CellContribution(a,t.cell_id,t.object)
  end
  for (id,t) in _mat.dict
    cellvals = get_cell_values(uhd,t.cell_id)
    a = attach_dirichlet(get_array(t),cellvals)
    matvec[id] = CellContribution(a,t.cell_id,t.object)
  end
  map(CollectionOfCellContribution, (matvec, mat, _vec.dict))
end


# AffineOperator

function AffineFEOperator(
  trial::FESpace,
  test::FESpace,
  assem::Assembler,
  biform::Function,
  liform::Function)

  u = get_cell_basis(trial)
  v = get_cell_basis(test)
  @assert is_trial(u)

  uhd = zero(trial)

  mat_terms = biform(u,v)
  vec_terms = liform(v)
  data = collect_cell_matrix_and_vector(mat_terms,vec_terms,uhd)
  A,b = assemble_matrix_and_vector(assem,data)

  AffineFEOperator(trial,test,A,b)
end

function AffineFEOperator(
  trial::FESpace,
  test::FESpace,
  biform::Function,
  liform::Function)

  assem = SparseMatrixAssembler(trial,test)
  AffineFEOperator(trial,test,assem,biform,liform)
end

function AffineFEOperator(
  mat::Type{<:AbstractSparseMatrix},
  trial::FESpace,
  test::FESpace,
  biform::Function,
  liform::Function)

  assem = SparseMatrixAssembler(mat,trial,test)
  AffineFEOperator(trial,test,assem,biform,liform)
end

function AffineFEOperator(
  mat::Type{<:AbstractSparseMatrix},
  vec::Type{<:AbstractVector},
  trial::FESpace,
  test::FESpace,
  biform::Function,
  liform::Function)

  assem = SparseMatrixAssembler(mat,vec,trial,test)
  AffineFEOperator(trial,test,assem,biform,liform)
end

# Some syntax sugar

struct BilinearForm
  a::Function
  U::FESpace
  V::FESpace
end

struct LinearForm
  b::Function
  V::FESpace
end

function Helpers.operate(a::Function,U::FESpace,V::FESpace)
  BilinearForm(a,U,V)
end

function Helpers.operate(b::Function,V::FESpace)
  LinearForm(b,V)
end

function Base.:(==)(bi::BilinearForm,li::LinearForm)
  @assert li.V === bi.V
  AffineFEOperator(bi.U,bi.V,bi.a,li.b)
end

function Base.:(==)(li::LinearForm,bi::BilinearForm)
  bi == li
end

"""
"""
macro form(fundef)
  s = "The @form macro is only allowed in function definitions"
  @assert isa(fundef,Expr) s
  @assert fundef.head in (:(=), :function) s
  funname = fundef.args[1].args[1]
  nargs = length(fundef.args[1].args)-1
  if nargs == 2

    U = fundef.args[1].args[2]
    V =  fundef.args[1].args[3]
    a = fundef.args[1].args[2:end]
    q = quote
      function $(funname)($(U)::FESpace,$(V)::FESpace)
        operate($(funname),$(a...))
      end
      $(fundef)
    end

  elseif nargs == 1

    V = fundef.args[1].args[2]
    a = fundef.args[1].args[2:end]
    q = quote
      function $(funname)($(V)::FESpace)
        operate($(funname),$(a...))
      end
      $(fundef)
    end

  else
    @unreachable "The @form macro is only allowed in functions with one or two arguments"
  end

  :($(esc(q)))
end
