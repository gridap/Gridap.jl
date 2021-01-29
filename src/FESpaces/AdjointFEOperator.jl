"""
AffineFEOperatorWithAdjoint

Reprepresent a fully assembled affine (linear) finite element problem
together with its corresponding adjoint operator.
See also [FEOperator](@ref)
"""
struct AffineFEOperatorWithAdjoint <: FEOperator
  trial::FESpace
  test::FESpace
  op::AffineOperatorWithAdjoint
end

"""
"""
function AffineFEOperatorWithAdjoint(
  trial::FESpace,
  test::FESpace,
  matrix::AbstractMatrix,
  matrix_adj::AbstractMatrix,
  vector::AbstractVector,
  vector_adj::AbstractVector)
  @assert num_free_dofs(trial) == size(matrix,2) "Incompatible trial space and matrix"
  @assert num_free_dofs(test) == size(matrix,1) "Incompatible test space and matrix"
  op = AffineOperatorWithAdjoint(matrix,matrix_adj,vector,vector_adj)
  AffineFEOperatorWithAdjoint(trial,test,op)
end

function AffineFEOperatorWithAdjoint(weakform::Function,assem::Assembler)

  trial = get_trial(assem)
  test = get_test(assem)

  u = get_cell_shapefuns_trial(trial)
  v = get_cell_shapefuns(test)

  uhd = zero(trial)
  matcontribs, matcontribs_adj, veccontribs, veccontribs_adj = weakform(u,v)
  data = collect_cell_matrix_and_vector(matcontribs,veccontribs,uhd)
  data_adj = collect_cell_matrix_and_vector(matcontribs_adj,veccontribs_adj)#,uhd) # Not sure about BCs...
  A,b = assemble_matrix_and_vector(assem,data)

  # @oriolcg: Not sure how to optimize this. The matrix A_adj doesn't need to be
  # assembled if we are using the LUSolver (it reuses the factorization of A).
  A_adj,b_adj = assemble_matrix_and_vector(assem,data_adj)

  AffineFEOperatorWithAdjoint(trial,test,A,A_adj,b,b_adj)
end

function AffineFEOperatorWithAdjoint(weakform::Function,args...)
  assem = SparseMatrixAssembler(args...)
  AffineFEOperatorWithAdjoint(weakform,assem)
end

function AffineFEOperatorWithAdjoint(a::Function,ℓ::Function,ℓd::Function,args...)
  AffineFEOperatorWithAdjoint(args...) do u,v
    a(u,v), a(v,u), ℓ(v), ℓd(u)
  end
end


get_test(feop::AffineFEOperatorWithAdjoint) = feop.test
get_trial(feop::AffineFEOperatorWithAdjoint) = feop.trial
get_matrix(feop::AffineFEOperatorWithAdjoint) = get_matrix(feop.op)
get_vector(feop::AffineFEOperatorWithAdjoint) = get_vector(feop.op)
get_adjoint_matrix(feop::AffineFEOperatorWithAdjoint) = get_adjoint_matrix(feop.op)
get_adjoint_vector(feop::AffineFEOperatorWithAdjoint) = get_adjoint_vector(feop.op)
get_algebraic_operator(feop::AffineFEOperatorWithAdjoint) = feop.op

function get_forward(feop::AffineFEOperatorWithAdjoint)
  test = get_test(feop)
  trial = get_trial(feop)
  matrix = get_matrix(feop)
  vector = get_vector(feop)
  op = AffineFEOperator(trial,test,matrix,vector)
end

function get_adjoint(feop::AffineFEOperatorWithAdjoint)
  test = get_trial(feop)
  trial = get_test(feop)
  matrix = get_adjoint_matrix(feop)
  vector = get_adjoint_vector(feop)
  op = AffineFEOperator(trial,test,matrix,vector)
end


# function AdjointFEOperator(op::FEOperator,xh::FEFunction)
#   @notimplemented
# end

# function AdjointFEOperatorWithAdjoint(op::AffineFEOperator,xh::FEFunction)
#   test = get_test(op)
#   trial = get_trial(op)
#   matrix = sparse(transpose(get_matrix(op)))
#   vector = allocate_residual(op,xh)
#   AffineFEOperator(test,trial,matrix,vector)
# end

# function AdjointFEOperator(op::FEOperatorFromWeakForm,xh::FEFunction)
#   test = get_test(op)
#   trial = get_trial(op)
#   dx = get_cell_shapefuns_trial(get_trial(op))
#   v = get_cell_shapefuns(get_test(op))
#   matdata = collect_cell_matrix(op.jac(xh,v,dx))
#   matrix = allocate_matrix(op.assem, matdata)
#   vector = allocate_residual(op,xh)
#   assemble_matrix!(matrix,op.assem,matdata)
#   AffineFEOperator(test,trial,matrix,vector)
# end
