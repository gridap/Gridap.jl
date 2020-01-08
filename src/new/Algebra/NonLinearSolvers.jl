
"""
    abstract type NonLinearSolver <: GridapType end

- [`solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator)`](@ref)
- [`solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator, cache)`](@ref)
"""
abstract type NonLinearSolver <: GridapType end

"""
    solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator)

Usage:

   cache = solve!(x,nls,op)

The returned `cache` object can be used in subsequent solves:

   solve!(x,nls,op,cache)
"""
function solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator)
  @abstractmethod
end

"""
    solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator,cache)

Solve using the cache object from a previous solve.
"""
function solve!(x::AbstractVector,nls::NonLinearSolver,op::NonLinearOperator,cache)
  @abstractmethod
end


"""
    solve(nls::NonLinearSolver,op::NonLinearOperator)

Creates and uses a zero initial guess.
"""
function solve(nls::NonLinearSolver,op::NonLinearOperator)
  x = allocate_solution(op)
  fill!(x,zero(eltype(x)))
  solve!(x,nls,op)
  x
end

"""
    test_non_linear_solver(
      nls::NonLinearSolver,
      op::NonLinearOperator,
      x0::AbstractVector,
      x::AbstractVector,
      pred::Function=isapprox)
"""
function test_non_linear_solver(
  nls::NonLinearSolver,
  op::NonLinearOperator,
  x0::AbstractVector,
  x::AbstractVector,
  pred::Function=isapprox)

  x1 = copy(x0)
  cache = solve!(x1,nls,op)
  @test pred(x1,x)

  x1 = copy(x0)
  solve!(x1,nls,op,cache)
  @test pred(x1,x)

end

# Mock

"""
    struct NewtonRaphsonSolver <:NonLinearSolver
      # Private fields
    end

Vanilla Newton-Raphson method
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



