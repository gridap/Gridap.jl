
struct FEOperatorFromWeakForm <: FEOperator
  res::Function
  jac::Function
  assem::Assembler
end

function FEOperator(res::Function,jac::Function,assem::Assembler)
  FEOperatorFromWeakForm(res,jac,assem)
end

function FEOperator(res::Function,jac::Function,args...)
  assem = SparseMatrixAssembler(args...)
  FEOperator(res,jac,assem)
end

function FEOperator(res::Function,assem::Assembler)
  jac(u,du,dv) = jacobian(x->res(x,dv),u)
  FEOperatorFromWeakForm(res,jac,assem)
end

function FEOperator(res::Function,args...)
  assem = SparseMatrixAssembler(args...)
  FEOperator(res,assem)
end

get_test(op::FEOperatorFromWeakForm) = get_test(op.assem)
get_trial(op::FEOperatorFromWeakForm) = get_trial(op.assem)

function allocate_residual(op::FEOperatorFromWeakForm,uh::FEFunction)
  v = get_cell_shapefuns(get_test(op))
  vecdata = collect_cell_vector(op.res(uh,v))
  allocate_vector(op.assem, vecdata)
end

function residual!(b::AbstractVector,op::FEOperatorFromWeakForm,uh::FEFunction)
  v = get_cell_shapefuns(get_test(op))
  vecdata = collect_cell_vector(op.res(uh,v))
  assemble_vector!(b,op.assem, vecdata)
  b
end

function allocate_jacobian(op::FEOperatorFromWeakForm,uh::FEFunction)
  du = get_cell_shapefuns_trial(get_trial(op))
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(op.jac(uh,du,v))
  allocate_matrix(op.assem, matdata)
end

function jacobian!(A::AbstractMatrix,op::FEOperatorFromWeakForm,uh::FEFunction)
  du = get_cell_shapefuns_trial(get_trial(op))
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(op.jac(uh,du,v))
  assemble_matrix!(A,op.assem,matdata)
  A
end

function residual_and_jacobian!(
  b::AbstractVector,A::AbstractMatrix,op::FEOperatorFromWeakForm,uh::FEFunction)
  du = get_cell_shapefuns_trial(get_trial(op))
  v = get_cell_shapefuns(get_test(op))
  data = collect_cell_matrix_and_vector(op.jac(uh,du,v),op.res(uh,v))
  assemble_matrix_and_vector!(A, b, op.assem, data)
  (b,A)
end

function residual_and_jacobian(op::FEOperatorFromWeakForm,uh::FEFunction)
  du = get_cell_shapefuns_trial(get_trial(op))
  v = get_cell_shapefuns(get_test(op))
  data = collect_cell_matrix_and_vector(op.jac(uh,du,v),op.res(uh,v))
  A, b = assemble_matrix_and_vector(op.assem, data)
  (b, A)
end

