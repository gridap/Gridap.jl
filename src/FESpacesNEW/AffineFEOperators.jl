
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

"""
   AffineFEOperator(test::FESpace,trial::FESpace,assem::Assembler,terms::AffineFETerm...)
   AffineFEOperator(test::FESpace,trial::FESpace,terms::AffineFETerm...)
"""
function AffineFEOperator(trial::FESpace,test::FESpace,assem::Assembler,terms::AffineFETerm...)

  u = get_cell_basis(trial)
  v = get_cell_basis(test)
  @assert is_trial(u)

  uhd = zero(trial)

  data = collect_cell_matrix_and_vector(uhd,u,v,terms)
  A,b = assemble_matrix_and_vector(assem,data)

  #matdata = collect_cell_matrix(u,v,terms)
  #vecdata = collect_cell_vector(uhd,v,terms)
  #A = assemble_matrix(assem,matdata)
  #b = assemble_vector(assem,vecdata)

  AffineFEOperator(trial,test,A,b)
end

function AffineFEOperator(trial::FESpace,test::FESpace,terms::AffineFETerm...)
  assem = SparseMatrixAssembler(trial,test)
  AffineFEOperator(trial,test,assem,terms...)
end

function AffineFEOperator(
  mat::Type{<:AbstractSparseMatrix},trial::FESpace,test::FESpace,terms::AffineFETerm...)
  assem = SparseMatrixAssembler(mat,trial,test)
  AffineFEOperator(trial,test,assem,terms...)
end

function AffineFEOperator(
  mat::Type{<:AbstractSparseMatrix},vec::Type{<:AbstractVector},trial::FESpace,test::FESpace,terms::AffineFETerm...)
  assem = SparseMatrixAssembler(mat,vec,trial,test)
  AffineFEOperator(trial,test,assem,terms...)
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

function allocate_residual(feop::AffineFEOperator,u)
  @assert is_a_fe_function(u)
  x = get_free_values(u)
  allocate_residual(feop.op,x)
end

function residual!(b::AbstractVector,feop::AffineFEOperator,u)
  @assert is_a_fe_function(u)
  x = get_free_values(u)
  residual!(b,feop.op,x)
end

function residual(feop::AffineFEOperator,u)
  @assert is_a_fe_function(u)
  x = get_free_values(u)
  residual(feop.op,x)
end

function allocate_jacobian(feop::AffineFEOperator,u)
  @assert is_a_fe_function(u)
  x = get_free_values(u)
  allocate_jacobian(feop.op,x)
end

function jacobian!(A::AbstractMatrix,feop::AffineFEOperator,u)
  @assert is_a_fe_function(u)
  x = get_free_values(u)
  jacobian!(A,feop.op,x)
end

function jacobian(feop::AffineFEOperator,u)
  @assert is_a_fe_function(u)
  x = get_free_values(u)
  jacobian(feop.op,x)
end
