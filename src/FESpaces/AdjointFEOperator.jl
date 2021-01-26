function AdjointFEOperator(op::FEOperator,xh::FEFunction)
  @notimplemented
end

function AdjointFEOperator(op::AffineFEOperator,xh::FEFunction)
  test = get_test(op)
  trial = get_trial(op)
  matrix = transpose(get_matrix(op))
  vector = zero(eltype(get_vector(op)))
  AffineFEOperator(test,trial,matrix,vector)
end

