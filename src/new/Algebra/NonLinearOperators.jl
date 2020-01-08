"""
    abstract type NonLinearOperator <: GridapType end

- [`residual!(b::AbstractVector,op::NonLinearOperator,x::AbstractVector)`](@ref)
- [`jacobian!(A::AbstractMatrix,op::NonLinearOperator,x::AbstractVector)`](@ref)
- [`num_domain_dims(op::NonLinearOperator)`](@ref)
- [`num_range_dims(op::NonLinearOperator)`](@ref)
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
    num_domain_dims(op::NonLinearOperator)
"""
function num_domain_dims(op::NonLinearOperator)
  @abstractmethod
end

"""
    num_range_dims(op::NonLinearOperator)
"""
function num_range_dims(op::NonLinearOperator)
  @abstractmethod
end

"""
    allocate_solution(op::NonLinearOperator)

"""
function allocate_solution(op::NonLinearOperator)
  Vector{Float64}(undef,num_domain_dims(op))
end

"""
    allocate_solution(op::NonLinearOperator,x::AbstractVector)

Allocate a solution with vector type equal to the given vector `x`.
"""
function allocate_solution(op::NonLinearOperator,x::AbstractVector)
  similar(x,eltype(x),num_domain_dims(op))
end

"""
    zero_initial_guess(op::NonLinearOperator)

Defaults to

    zeros(num_domain_dims(op))
"""
function zero_initial_guess(op::NonLinearOperator)
  zeros(num_domain_dims(op))
end

"""
    zero_initial_guess(op::NonLinearOperator,x::AbstractVector)

The default implementation of this function is

    similar(x,eltype(x),num_domain_dims(op))
    fill!(x,zero(eltype(x)))
"""
function zero_initial_guess(op::NonLinearOperator,x::AbstractVector)
  similar(x,eltype(x),num_domain_dims(op))
  fill!(x,zero(eltype(x)))
  x
end

"""
    allocate_residual(op::NonLinearOperator,x::AbstractVector)
"""
function allocate_residual(op::NonLinearOperator,x::AbstractVector)
  similar(x,eltype(x),num_range_dims(op))
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

  @test length(b) == num_range_dims(op)
  @test length(x) == num_domain_dims(op)
  b1 = allocate_residual(op,x)
  residual!(b1,op,x)
  @test pred(b,b1)

  if jac != nothing
    nrows, ncols = size(jac)
    @test nrows == num_range_dims(op)
    @test ncols == num_domain_dims(op)
    A = allocate_jacobian(op,x)
    jacobian!(A,op,x)
    @test pred(A,jac)
  end

end

"""
    is_square(op::NonLinearOperator)
"""
function is_square(op::NonLinearOperator)
  num_range_dims(op) == num_domain_dims(op)
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

num_range_dims(op::NonLinearOperatorMock) = 2

num_domain_dims(op::NonLinearOperatorMock) = num_range_dims(op)

function allocate_jacobian(op::NonLinearOperatorMock,x::AbstractVector)
  T = eltype(x)
  n = num_range_dims(op)
  m = num_domain_dims(op)
  zeros(T,n,m)
end

