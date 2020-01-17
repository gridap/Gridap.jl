
"""
"""
abstract type FEOperator <: NonLinearOperator end

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

function allocate_residual(op::FEOperator,x::AbstractVector)
  trial = get_trial(op)
  u = FEFunction(trial,x)
  allocate_residual(op,u)
end

"""
"""
function residual!(b::AbstractVector,op::FEOperator,u)
  @assert is_a_fe_function(u)
  @abstractmethod
end

function residual!(b::AbstractVector,op::FEOperator,x::AbstractVector)
  trial = get_trial(op)
  u = FEFunction(trial,x)
  residual!(b,op,u)
end

function residual(op::FEOperator,u)
  @assert is_a_fe_function(u)
  b = allocate_residual(op,u)
  residual!(b,op,u)
  b
end

function residual(op::FEOperator,x::AbstractVector)
  trial = get_trial(op)
  u = FEFunction(trial,x)
  residual(op,u)
end

"""
"""
function allocate_jacobian(op::FEOperator,u)
  @assert is_a_fe_function(u)
  @abstractmethod
end

function allocate_jacobian(op::FEOperator,x::AbstractVector)
  trial = get_trial(op)
  u = FEFunction(trial,x)
  allocate_jacobian(op,u)
end

"""
"""
function jacobian!(A::AbstractMatrix,op::FEOperator,u)
  @assert is_a_fe_function(u)
  @abstractmethod
end

function jacobian!(A::AbstractMatrix,op::FEOperator,x::AbstractVector)
  trial = get_trial(op)
  u = FEFunction(trial,x)
  jacobian!(A,op,u)
end

function jacobian(op::FEOperator,u)
  @assert is_a_fe_function(u)
  A = allocate_jacobian(op,u)
  jacobian!(A,op,u)
  A
end

function jacobian(op::FEOperator,x::AbstractVector)
  trial = get_trial(op)
  u = FEFunction(trial,x)
  jacobian(op,u)
end

function zero_initial_guess(::Type{T},op::FEOperator) where T
  trial = get_trial(op)
  x = zero_free_values(T,trial)
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
  test_non_linear_operator(op,args...;kwargs...)
end

