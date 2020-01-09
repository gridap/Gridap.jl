"""
    abstract type NonLinearOperator <: GridapType end

- [`residual!(b::AbstractVector,op::NonLinearOperator,x::AbstractVector)`](@ref)
- [`jacobian!(A::AbstractMatrix,op::NonLinearOperator,x::AbstractVector)`](@ref)
- [`zero_initial_guess(op::NonLinearOperator)`](@ref)
- [`zero_initial_guess(op::NonLinearOperator,x::AbstractVector)`](@ref)
- [`allocate_residual(op::NonLinearOperator,x::AbstractVector)`](@ref)
- [`allocate_jacobian(op::NonLinearOperator,x::AbstractVector)`](@ref)

"""
abstract type NonLinearOperator <: GridapType end

"""
    residual!(b::AbstractVector,op::NonLinearOperator,x::AbstractVector)
"""
function residual!(b::AbstractVector,op::NonLinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    residual(op::NonLinearOperator,x::AbstractVector)
"""
function residual(op::NonLinearOperator,x::AbstractVector)
  b = allocate_residual(op,x)
  residual!(b,op,x)
  b
end

"""
    jacobian!(A::AbstractMatrix,op::NonLinearOperator,x::AbstractVector)
"""
function jacobian!(A::AbstractMatrix,op::NonLinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    jacobian(op::NonLinearOperator,x::AbstractVector)
"""
function jacobian(op::NonLinearOperator,x::AbstractVector)
  A = allocate_jacobian(op,x)
  jacobian!(A,op,x)
  A
end

"""
    residual_and_jacobian!(
      b::AbstractVector, A::AbstractMatrix,
      op::NonLinearOperator, x::AbstractVector)
"""
function residual_and_jacobian!(
  b::AbstractVector, A::AbstractMatrix, op::NonLinearOperator, x::AbstractVector)
  residual!(b,op,x)
  jacobian!(A,op,x)
  (b, A)
end

"""
    residual_and_jacobian(op::NonLinearOperator,x::AbstractVector)
"""
function residual_and_jacobian(op::NonLinearOperator,x::AbstractVector)
  b, A = allocate_residual_and_jacobian(op,x)
  residual_and_jacobian!(b,A,op,x)
  (b, A)
end

"""
    zero_initial_guess(op::NonLinearOperator)
"""
function zero_initial_guess(op::NonLinearOperator)
  @abstractmethod
end

"""
    zero_initial_guess(op::NonLinearOperator,x::AbstractVector)
"""
function zero_initial_guess(op::NonLinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    allocate_residual(op::NonLinearOperator,x::AbstractVector)
"""
function allocate_residual(op::NonLinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    allocate_jacobian(op::NonLinearOperator,x::AbstractVector)
"""
function allocate_jacobian(op::NonLinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    allocate_residual_and_jacobian(op::NonLinearOperator,x::AbstractVector)
"""
function allocate_residual_and_jacobian(op::NonLinearOperator,x::AbstractVector)
  b = allocate_residual(op,x)
  A = allocate_jacobian(op,x)
  (b,A)
end

"""
    test_non_linear_operator(
      op::NonLinearOperator,
      x::AbstractVector,
      b::AbstractVector,
      pred=isapprox;
      jac=nothing)
"""
function test_non_linear_operator(
  op::NonLinearOperator,
  x::AbstractVector,
  b::AbstractVector,
  pred=isapprox;
  jac=nothing)

  b1 = allocate_residual(op,x)
  residual!(b1,op,x)
  @test pred(b,b1)

  x0 = zero_initial_guess(op)
  x0 = zero_initial_guess(op,x)

  if jac != nothing
    nrows, ncols = size(jac)
    A = allocate_jacobian(op,x)
    jacobian!(A,op,x)
    @test pred(A,jac)
  end

end

# Mock for testing purposes

struct NonLinearOperatorMock <: NonLinearOperator end

function residual!(b::AbstractVector,::NonLinearOperatorMock,x::AbstractVector)
  b[1] = (x[1]-2)*(x[1]-1)
  b[2] = x[2]-3
end

function jacobian!(A::AbstractMatrix,::NonLinearOperatorMock,x::AbstractVector)
  A[1,1] = (x[1]-1) + (x[1]-2)
  A[1,2] = 0
  A[2,1] = 0
  A[2,2] = 1
end

function zero_initial_guess(op::NonLinearOperatorMock)
  x = Float64[]
  allocate_residual(op,x)
end

function zero_initial_guess(op::NonLinearOperatorMock,x::AbstractVector)
  allocate_residual(op,x)
end

function allocate_residual(op::NonLinearOperatorMock,x::AbstractVector)
  T = eltype(x)
  n = 2
  zeros(T,n)
end

function allocate_jacobian(op::NonLinearOperatorMock,x::AbstractVector)
  T = eltype(x)
  n = 2
  m = 2
  zeros(T,n,m)
end

