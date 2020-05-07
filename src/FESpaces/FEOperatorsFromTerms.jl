
struct FEOperatorFromTerms <: FEOperator
  trial::FESpace
  test::FESpace
  assem::Assembler
  terms
end

"""
"""
function FEOperator(trial::FESpace,test::FESpace,assem::Assembler,terms::FETerm...)
  FEOperatorFromTerms(trial,test,assem,terms)
end

"""
"""
function FEOperator(trial::FESpace,test::FESpace,terms::FETerm...)
  assem = SparseMatrixAssembler(trial,test)
  FEOperator(trial,test,assem,terms...)
end

function FEOperator(mat::Type{<:AbstractSparseMatrix},trial::FESpace,test::FESpace,terms::FETerm...)
  assem = SparseMatrixAssembler(mat,trial,test)
  FEOperator(trial,test,assem,terms...)
end

function FEOperator(mat::Type{<:AbstractSparseMatrix},vec::Type{<:AbstractVector},trial::FESpace,test::FESpace,terms::FETerm...)
  assem = SparseMatrixAssembler(mat,vec,trial,test)
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
  vecdata = collect_cell_residual(uh,v,op.terms)
  allocate_vector(op.assem, vecdata)
end

function residual!(b::AbstractVector,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  v = get_cell_basis(op.test)
  vecdata = collect_cell_residual(uh,v,op.terms)
  assemble_vector!(b,op.assem, vecdata)
  b
end

function allocate_jacobian(op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  matdata = collect_cell_jacobian(uh,du,v,op.terms)
  allocate_matrix(op.assem, matdata)
end

function jacobian!(A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  matdata = collect_cell_jacobian(uh,du,v,op.terms)
  assemble_matrix!(A,op.assem,matdata)
  A
end

function residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  data = collect_cell_jacobian_and_residual(uh,du,v,op.terms)
  assemble_matrix_and_vector!(A, b, op.assem, data)
  (b,A)
end

function residual_and_jacobian(op::FEOperatorFromTerms,uh)
  @assert is_a_fe_function(uh)
  du = get_cell_basis(op.trial)
  v = get_cell_basis(op.test)
  data = collect_cell_jacobian_and_residual(uh,du,v,op.terms)
  A, b = assemble_matrix_and_vector(op.assem, data)
  (b, A)
end
