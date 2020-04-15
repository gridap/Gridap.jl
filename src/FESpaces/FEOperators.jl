
"""
"""
abstract type FEOperator <: GridapType end

"""
"""
function get_test(op::FEOperator)
  @abstractmethod
end

"""
"""
function get_trial(op::FEOperator)
  @abstractmethod
end

"""
"""
function allocate_residual(op::FEOperator,u)
  @assert is_a_fe_function(u)
  @abstractmethod
end

"""
"""
function residual!(b::AbstractVector,op::FEOperator,u)
  @assert is_a_fe_function(u)
  @abstractmethod
end

function residual(op::FEOperator,u)
  @assert is_a_fe_function(u)
  b = allocate_residual(op,u)
  residual!(b,op,u)
  b
end

"""
"""
function allocate_jacobian(op::FEOperator,u)
  @assert is_a_fe_function(u)
  @abstractmethod
end

"""
"""
function jacobian!(A::AbstractMatrix,op::FEOperator,u)
  @assert is_a_fe_function(u)
  @abstractmethod
end

function jacobian(op::FEOperator,u)
  @assert is_a_fe_function(u)
  A = allocate_jacobian(op,u)
  jacobian!(A,op,u)
  A
end

"""
"""
function residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::FEOperator,u)
  residual!(b,op,u)
  jacobian!(A,op,u)
  (b,A)
end

"""
"""
function residual_and_jacobian(op::FEOperator,u)
  b = residual(op,u)
  A = jacobian(op,u)
  (b,A)
end

"""
"""
function test_fe_operator(op::FEOperator,args...;kwargs...)
  test = get_test(op)
  trial = get_trial(op)
  @test isa(test,FESpace)
  @test isa(trial,FESpace)
  u = zero(trial)
  @test is_a_fe_function(u)
  b = allocate_residual(op,u)
  @test isa(b,AbstractVector)
  residual!(b,op,u)
  b2 = residual(op,u)
  @test isa(b2,AbstractVector)
  A = allocate_jacobian(op,u)
  @test isa(A,AbstractMatrix)
  jacobian!(A,op,u)
  A2 = jacobian(op,u)
  @test isa(A2,AbstractMatrix)
  residual_and_jacobian!(b,A,op,u)
  b, A = residual_and_jacobian(op,u)
  @test isa(b,AbstractVector)
  @test isa(A,AbstractMatrix)
  _op = get_algebraic_operator(op)
  test_nonlinear_operator(_op,args...;kwargs...)
end

# FEOperator viewed as a NonlinearOperator

"""
"""
function get_algebraic_operator(feop::FEOperator)
  AlgebraicOpFromFEOp(feop)
end

struct AlgebraicOpFromFEOp <: NonlinearOperator
  feop::FEOperator
end

function allocate_residual(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  allocate_residual(op.feop,u)
end

function residual!(b::AbstractVector,op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual!(b,op.feop,u)
end

function residual(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual(op.feop,u)
end

function allocate_jacobian(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  allocate_jacobian(op.feop,u)
end

function jacobian!(A::AbstractMatrix,op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  jacobian!(A,op.feop,u)
end

function jacobian(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  jacobian(op.feop,u)
end

function residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual_and_jacobian!(b,A,op.feop,u)
end

function residual_and_jacobian(op::AlgebraicOpFromFEOp,x::AbstractVector)
  trial = get_trial(op.feop)
  u = EvaluationFunction(trial,x)
  residual_and_jacobian(op.feop,u)
end

function zero_initial_guess(::Type{T},op::AlgebraicOpFromFEOp) where T
  trial = get_trial(op.feop)
  x = zero_free_values(T,trial)
end
