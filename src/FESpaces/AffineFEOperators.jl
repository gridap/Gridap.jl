
"""
    AffineFEOperator

Reprepresent a fully assembled affine (linear) finite element problem.
See also [FEOperator](@ref)
"""
struct AffineFEOperator <: FEOperator
  trial::FESpace
  test::FESpace
  op::AffineOperator
end

"""
"""
function AffineFEOperator(trial::FESpace,test::FESpace,matrix::AbstractMatrix,vector::AbstractVector)
  @assert num_free_dofs(trial) == size(matrix,2) "Incompatible trial space and matrix"
  @assert num_free_dofs(test) == size(matrix,1) "Incompatible test space and matrix"
  op = AffineOperator(matrix,vector)
  AffineFEOperator(trial,test,op)
end

function AffineFEOperator(weakform::Function,assem::Assembler)

  trial = get_trial(assem)
  test = get_test(assem)

  u = get_cell_shapefuns_trial(trial)
  v = get_cell_shapefuns(test)

  uhd = zero(trial)
  matcontribs, veccontribs = weakform(u,v)
  data = collect_cell_matrix_and_vector(matcontribs,veccontribs,uhd)
  A,b = assemble_matrix_and_vector(assem,data)

  AffineFEOperator(trial,test,A,b)
end

function AffineFEOperator(weakform::Function,args...)
  assem = SparseMatrixAssembler(args...)
  AffineFEOperator(weakform,assem)
end

function AffineFEOperator(a::Function,ℓ::Function,args...)
  AffineFEOperator(args...) do u,v
    a(u,v), ℓ(v)
  end
end

get_test(feop::AffineFEOperator) = feop.test

get_trial(feop::AffineFEOperator) = feop.trial

"""
"""
get_matrix(feop::AffineFEOperator) = get_matrix(feop.op)

"""
"""
get_vector(feop::AffineFEOperator) = get_vector(feop.op)

get_algebraic_operator(feop::AffineFEOperator) = feop.op

function allocate_residual(feop::AffineFEOperator,u::FEFunction)
  x = get_free_values(u)
  allocate_residual(feop.op,x)
end

function residual!(b::AbstractVector,feop::AffineFEOperator,u::FEFunction)
  x = get_free_values(u)
  residual!(b,feop.op,x)
end

function residual(feop::AffineFEOperator,u::FEFunction)
  x = get_free_values(u)
  residual(feop.op,x)
end

function allocate_jacobian(feop::AffineFEOperator,u::FEFunction)
  x = get_free_values(u)
  allocate_jacobian(feop.op,x)
end

function jacobian!(A::AbstractMatrix,feop::AffineFEOperator,u::FEFunction)
  x = get_free_values(u)
  jacobian!(A,feop.op,x)
end

function jacobian(feop::AffineFEOperator,u::FEFunction)
  x = get_free_values(u)
  jacobian(feop.op,x)
end
