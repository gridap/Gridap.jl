"""
    abstract type NonlinearOperator <: GridapType end

- [`residual!(b::AbstractVector,op::NonlinearOperator,x::AbstractVector)`](@ref)
- [`jacobian!(A::AbstractMatrix,op::NonlinearOperator,x::AbstractVector)`](@ref)
- [`zero_initial_guess(op::NonlinearOperator)`](@ref)
- [`allocate_residual(op::NonlinearOperator,x::AbstractVector)`](@ref)
- [`allocate_jacobian(op::NonlinearOperator,x::AbstractVector)`](@ref)

"""
abstract type NonlinearOperator <: GridapType end

"""
    residual!(b::AbstractVector,op::NonlinearOperator,x::AbstractVector)
"""
function residual!(b::AbstractVector,op::NonlinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    residual(op::NonlinearOperator,x::AbstractVector)
"""
function residual(op::NonlinearOperator,x::AbstractVector)
  b = allocate_residual(op,x)
  residual!(b,op,x)
  b
end

"""
    jacobian!(A::AbstractMatrix,op::NonlinearOperator,x::AbstractVector)
"""
function jacobian!(A::AbstractMatrix,op::NonlinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    jacobian(op::NonlinearOperator,x::AbstractVector)
"""
function jacobian(op::NonlinearOperator,x::AbstractVector)
  A = allocate_jacobian(op,x)
  jacobian!(A,op,x)
  A
end

"""
    residual_and_jacobian!(
      b::AbstractVector, A::AbstractMatrix,
      op::NonlinearOperator, x::AbstractVector)
"""
function residual_and_jacobian!(
  b::AbstractVector, A::AbstractMatrix, op::NonlinearOperator, x::AbstractVector)
  residual!(b,op,x)
  jacobian!(A,op,x)
  (b, A)
end

"""
    residual_and_jacobian(op::NonlinearOperator,x::AbstractVector)
"""
function residual_and_jacobian(op::NonlinearOperator,x::AbstractVector)
  b, A = allocate_residual_and_jacobian(op,x)
  residual_and_jacobian!(b,A,op,x)
  (b, A)
end

function hessian end
function hessian! end

"""
    zero_initial_guess(op::NonlinearOperator)
"""
function zero_initial_guess(op::NonlinearOperator)
  @abstractmethod
end

"""
    allocate_residual(op::NonlinearOperator,x::AbstractVector)
"""
function allocate_residual(op::NonlinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    allocate_jacobian(op::NonlinearOperator,x::AbstractVector)
"""
function allocate_jacobian(op::NonlinearOperator,x::AbstractVector)
  @abstractmethod
end

"""
    allocate_residual_and_jacobian(op::NonlinearOperator,x::AbstractVector)
"""
function allocate_residual_and_jacobian(op::NonlinearOperator,x::AbstractVector)
  b = allocate_residual(op,x)
  A = allocate_jacobian(op,x)
  (b,A)
end

"""
    test_nonlinear_operator(
      op::NonlinearOperator,
      x::AbstractVector,
      b::AbstractVector,
      pred=isapprox;
      jac=nothing)
"""
function test_nonlinear_operator(
  op::NonlinearOperator,
  x::AbstractVector,
  b::AbstractVector,
  pred=isapprox;
  jac=nothing)

  b1 = allocate_residual(op,x)
  residual!(b1,op,x)
  @test pred(b,b1)

  x0 = zero_initial_guess(op)

  if jac != nothing
    nrows, ncols = size(jac)
    A = allocate_jacobian(op,x)
    jacobian!(A,op,x)
    @test pred(A,jac)

    residual_and_jacobian!(b1,A,op,x)
    @test pred(A,jac)
    @test pred(b,b1)

    b1, A = residual_and_jacobian(op,x)
    @test pred(A,jac)
    @test pred(b,b1)
  end

end

# Mock for testing purposes

struct NonlinearOperatorMock <: NonlinearOperator end

function residual!(b::AbstractVector,::NonlinearOperatorMock,x::AbstractVector)
  b[1] = (x[1]-2)*(x[1]-1)
  b[2] = x[2]-3
  b
end

function jacobian!(A::AbstractMatrix,::NonlinearOperatorMock,x::AbstractVector)
  A[1,1] = (x[1]-1) + (x[1]-2)
  A[1,2] = 0
  A[2,1] = 0
  A[2,2] = 1
  A
end

function zero_initial_guess(op::NonlinearOperatorMock)
  x = allocate_residual(op,Float64[])
  fill!(x,zero(eltype(x)))
  x
end

function allocate_residual(op::NonlinearOperatorMock,x::AbstractVector)
  T = eltype(x)
  n = 2
  similar(x,T,n)
end

function allocate_jacobian(op::NonlinearOperatorMock,x::AbstractVector)
  T = eltype(x)
  n = 2
  m = 2
  zeros(T,n,m)
end
