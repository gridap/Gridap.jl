module NonLinearSolvers

using Gridap.LinearSolvers
using Gridap.Helpers

export NonLinearOperator
export residual
export jacobian
export residual!
export jacobian!
export create_in_domain

export NonLinearSolver
export NextNonLinearSolver
import Gridap: solve
import Gridap: solve!

export NewtonRaphsonSolver

abstract type NonLinearOperator end

function residual!(b::AbstractVector,::NonLinearOperator,x::AbstractVector)
  @abstractmethod
end

function jacobian!(A::AbstractMatrix,::NonLinearOperator,x::AbstractVector)
  @abstractmethod
end

function jacobian(::NonLinearOperator,x::AbstractVector)::AbstractMatrix
  @abstractmethod
end

function create_in_domain(::NonLinearOperator)::AbstractVector
  @abstractmethod
end

function residual(op::NonLinearOperator,x::AbstractVector)
  b = similar(x)
  residual!(b,op,x)
  b
end

abstract type NonLinearSolver end

function solve!(x::AbstractVector,::NonLinearSolver,::NonLinearOperator)::Any
  @abstractmethod
end

function solve!(x::AbstractVector,::NonLinearSolver,::NonLinearOperator,cache::Any)
  @abstractmethod
end

function solve(nls::NonLinearSolver,op::NonLinearOperator)
  x = create_in_domain(op)
  E = eltype(x)
  x .= zero(E)
  solve!(x,nls,op)
  x
end

"""
Newton-Raphson method
"""
struct NewtonRaphsonSolver <:NonLinearSolver
  ls::LinearSolver
  tol::Float64
  max_nliters::Int
end

struct NewtonRaphsonCache
  A::AbstractMatrix
  b::AbstractVector
  dx::AbstractVector
  ns::NumericalSetup
end

function solve!(x::AbstractVector,nls::NewtonRaphsonSolver,op::NonLinearOperator)

  b = residual(op, x)
  A = jacobian(op, x)
  dx = similar(b)
  ss = symbolic_setup(nls.ls, A)
  ns = numerical_setup(ss,A)

  _solve_nr!(x,A,b,dx,ns,nls,op)

  NewtonRaphsonCache(A,b,dx,ns)

end

function solve!(
  x::AbstractVector,nls::NewtonRaphsonSolver,op::NonLinearOperator,cache::NewtonRaphsonCache)

  b = cache.b
  A = cache.A
  dx = cache.dx
  ns = cache.ns

  residual!(b, op, x)
  jacobian!(A, op, x)
  numerical_setup!(ns,A)

  _solve_nr!(x,A,b,dx,ns,nls,op)

end

function _solve_nr!(x,A,b,dx,ns,nls,op)

  # Check for convergence on the initial residual
  isconv, conv0 = _check_convergence(nls,b)
  if isconv; return; end

  # Newton-like iterations
  for nliter in 1:nls.max_nliters

    # Solve linearized problem
    broadcast!(*,b,b,-1)
    solve!(dx,ns,b)
    broadcast!(+,x,x,dx)

    # Check convergence for the current residual
    residual!(b, op, x)
    isconv = _check_convergence(nls, b, conv0)
    if isconv; return; end

    if nliter == nls.max_nliters
      @unreachable
    end

    # Assemble jacobian (fast in-place version)
    # and prepare solver
    jacobian!(A, op, x)
    numerical_setup!(ns,A)

  end

end

function _check_convergence(nls,b)
  m0 = _inf_norm(b)
  (false, m0)
end

function _check_convergence(nls,b,m0)
  m = _inf_norm(b)
  m < nls.tol * m0
end

function _inf_norm(b)
  m = 0
  for bi in b
    m = max(m,abs(bi))
  end
  m
end

end # module NonLinearSolvers
