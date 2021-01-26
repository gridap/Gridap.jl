function AdjointFEOperator(op::FEOperator,xh::FEFunction)
  @notimplemented
end

function AdjointFEOperator(op::AffineFEOperator,xh::FEFunction)
  test = get_test(op)
  trial = get_trial(op)
  matrix = sparse(transpose(get_matrix(op)))
  vector = allocate_residual(op,xh)
  AffineFEOperator(test,trial,matrix,vector)
end

function AdjointFEOperator(op::FEOperatorFromWeakForm,xh::FEFunction)
  test = get_test(op)
  trial = get_trial(op)
  dx = get_cell_shapefuns_trial(get_trial(op))
  v = get_cell_shapefuns(get_test(op))
  matdata = collect_cell_matrix(op.jac(xh,v,dx))
  matrix = allocate_matrix(op.assem, matdata)
  vector = allocate_residual(op,xh)
  assemble_matrix!(matrix,op.assem,matdata)
  AffineFEOperator(test,trial,matrix,vector)
end
