
struct FEOperatorFromTerms <: FEOperator
  trial::FESpace
  test::FESpace
  assem::Assembler
  terms
  function FEOperatorFromTerms(trial::FESpace,test::FESpace,assem::Assembler,terms::FETerm...)
    new(trial,test,assem,terms)
  end
end

"""
"""
function FEOperator(trial::FESpace,test::FESpace,assem::Assembler,terms::FETerm...)
  FEOperatorFromTerms(trial,test,assem,terms...)
end

"""
"""
function FEOperator(trial::FESpace,test::FESpace,terms::FETerm...)
  assem = SparseMatrixAssembler(test,trial)
  FEOperator(trial,test,assem,terms...)
end

function FEOperator(mat::Type{<:AbstractSparseMatrix},trial::FESpace,test::FESpace,terms::FETerm...)
  assem = SparseMatrixAssembler(mat,test,trial)
  FEOperator(trial,test,assem,terms...)
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
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  _, cellidsrows, cellidscols = collect_cell_jacobian(uh,v,du,op.terms)
  allocate_matrix(op.assem, cellidsrows, cellidscols)
end

function jacobian!(A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  cellmats, cellidsrows, cellidscols = collect_cell_jacobian(uh,v,du,op.terms)
  assemble_matrix!(A,op.assem, cellmats, cellidsrows, cellidscols)
  A
end

function residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  data = collect_cell_jacobian_and_residual(uh,v,du,op.terms)
  assemble_matrix_and_vector!(A, b, op.assem,data...)
  (b,A)
end

function residual_and_jacobian(op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  data = collect_cell_jacobian_and_residual(uh,v,du,op.terms)
  A, b = assemble_matrix_and_vector(op.assem,data...)
  (b, A)
end
