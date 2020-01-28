
struct FEOperatorFromTerms <: FEOperator
  test::FESpace
  trial::FESpace
  assem::Assembler
  terms
  function FEOperatorFromTerms(test::FESpace,trial::FESpace,assem::Assembler,terms::FETerm...)
    new(test,trial,assem,terms)
  end
end

"""
"""
function FEOperator(test::FESpace,trial::FESpace,assem::Assembler,terms::FETerm...)
  FEOperatorFromTerms(test,trial,assem,terms...)
end

"""
"""
function FEOperator(test::FESpace,trial::FESpace,terms::FETerm...)
  assem = SparseMatrixAssembler(test,trial)
  FEOperator(test,trial,assem,terms...)
end

function FEOperator(mat::Type{<:AbstractSparseMatrix},test::FESpace,trial::FESpace,terms::FETerm...)
  assem = SparseMatrixAssembler(mat,test,trial)
  FEOperator(test,trial,assem,terms...)
end

function get_test(op::FEOperatorFromTerms)
  op.test
end

function get_trial(op::FEOperatorFromTerms)
  op.trial
end

function allocate_residual(op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  _, cellids = collect_cell_residual(uh,v,op.terms)
  allocate_vector(op.assem, cellids)
end

function residual!(b::AbstractVector,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  cellvecs, cellids = collect_cell_residual(uh,v,op.terms)
  assemble_vector!(b,op.assem, cellvecs, cellids)
  b
end

function allocate_jacobian(op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  du = get_cell_basis(op.trial)
  _, cellidsrows, cellidscols = collect_cell_jacobian(uh,v,du,op.terms)
  allocate_matrix(op.assem, cellidsrows, cellidscols)
end

function jacobian!(A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  du = get_cell_basis(op.trial)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian(uh,v,du,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

