module FEOperators

using Gridap.Helpers
using Gridap.Geometry
using Gridap.CellValues
using Gridap.CellMaps
using Gridap.CellQuadratures
using Gridap.CellIntegration

using Gridap.FESpaces
using Gridap.Assemblers
using Gridap.LinearSolvers
using LinearAlgebra

export FEOperator
export FESolver
export LinearFEOperator
export LinearFESolver
export NonLinearFEOperator
export NonLinearFESolver
import Gridap: solve
import Gridap: solve!
export jacobian
export jacobian!
import Gridap: apply
import Gridap: apply!
import Gridap.FESpaces: TrialFESpace
import Gridap.FESpaces: TestFESpace
import Gridap.LinearSolvers: LinearSolver

abstract type FEOperator end

function apply(::FEOperator,::FEFunction)::AbstractVector
  @abstractmethod
end

function jacobian(::FEOperator,::FEFunction)::AbstractMatrix
  @abstractmethod
end

function apply!(::AbstractVector,::FEOperator,::FEFunction)
  @abstractmethod
end

function jacobian!(::AbstractMatrix,::FEOperator,::FEFunction)
  @abstractmethod
end

function TrialFESpace(::FEOperator)::FESpaceWithDirichletData
  @abstractmethod
end

function TestFESpace(::FEOperator)::FESpaceWithDirichletData
  @abstractmethod
end

abstract type FESolver end

LinearSolver(::FESolver)::LinearSolver = @abstractmethod

"""
Solve with a given initial guess that will be overwritten with the solution.
The initial Jacobian, initial residual, and the initial numerical setup are also given.
This routine is useful is useful for time-dependent problems, where the initial
data has been already set in another place for all time steps.
"""
function solve!(
  uh::FEFunction,A::AbstractMatrix,b::AbstractVector,dx::AbstractVector,
  ::NumericalSetup,::FESolver,::FEOperator)
  @abstractmethod
end

"""
Solve with given initial guess, which is overwritten by the solution.
In contrast to previous one, this method does the set up of the first linear solve.
This routine is useful is for steady-state problems.
"""
function solve!(uh::FEFunction,nls::FESolver,op::FEOperator)

  # Setup system for first solve
  b = apply(op, uh)
  A = jacobian(op, uh)
  dx = similar(b)
  ls = LinearSolver(nls)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)

  # Do the iterations
  solve!(uh,A,b,dx,ns,nls,op)

end

"""
Solve that allocates, and sets initial guess to zero
and returns the solution
"""
function solve(nls::FESolver,op::FEOperator)

  # Setup initial condition
  U = TrialFESpace(op)
  uh = zero(U)

  # Solve
  solve!(uh,nls,op)
  uh
end

"""
Struct representing a linear FE Operator
"""
struct LinearFEOperator{M,V} <:FEOperator
  mat::M
  vec::V
  trialfesp::FESpaceWithDirichletData
  testfesp::FESpaceWithDirichletData
end

function LinearFEOperator(
  biform::Function,
  liform::Function,
  testfesp::FESpace{D,Z,T},
  trialfesp::FESpaceWithDirichletData,
  assem::Assembler{M,V},
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where {M,V,D,Z,T}

  # This will not be a CellBasis in the future
  v = CellBasis(testfesp)
  u = CellBasis(trialfesp)

  # The way we modify the rhs can be improved
  E = eltype(T)
  free_values = zeros(E,num_free_dofs(trialfesp))
  diri_values = diri_dofs(trialfesp)
  uhd = FEFunction(trialfesp,free_values,diri_values)

  cellmat = integrate(biform(v,u),trian,quad)
  cellvec = integrate( liform(v)-biform(v,uhd), trian, quad)

  mat = assemble(assem,cellmat)
  vec = assemble(assem,cellvec)

  LinearFEOperator(mat,vec,trialfesp,testfesp)

end

TrialFESpace(op::LinearFEOperator) = op.trialfesp

TestFESpace(op::LinearFEOperator) = op.testfesp

function apply(o::LinearFEOperator,uh::FEFunction)
  vals = free_dofs(uh)
  o.mat * vals - o.vec
end

function apply!(b::Vector,o::LinearFEOperator,uh::FEFunction)
  vals = free_dofs(uh)
  mul!(b,o.mat,vals)
  broadcast!(-,b,b,o.vec)
end

jacobian(o::LinearFEOperator,::FEFunction) = o.mat

function jacobian!(mat::AbstractMatrix, o::LinearFEOperator, ::FEFunction)
  @assert mat === o.mat
end

"""
The solver that solves a LinearFEOperator
"""
struct LinearFESolver <: FESolver
  ls::LinearSolver
end

LinearSolver(s::LinearFESolver) = s.ls

function solve!(
  ::FEFunction,::AbstractMatrix,::AbstractVector,
  ::NumericalSetup,::LinearFESolver,::FEOperator)
  @unreachable
end

function solve!(
  uh::FEFunction,A::AbstractMatrix,b::AbstractVector,dx::AbstractVector,ns::NumericalSetup,
  s::LinearFESolver,o::LinearFEOperator)

  x = free_dofs(uh)
  broadcast!(*,b,b,-1)
  solve!(dx,ns,A,b)
  broadcast!(+,x,x,dx)
end

"""
Struct representing a non-linear FE Operator
"""
struct NonLinearFEOperator{D,Z,T} <:FEOperator
  res::Function
  jac::Function
  testfesp::FESpaceWithDirichletData{D,Z,T}
  trialfesp::FESpaceWithDirichletData{D,Z,T}
  assem::Assembler
  trian::Triangulation{Z}
  quad::CellQuadrature{Z}
end

TrialFESpace(op::NonLinearFEOperator) = op.trialfesp

TestFESpace(op::NonLinearFEOperator) = op.testfesp

function apply(op::NonLinearFEOperator,uh::FEFunction)
  cellvec = _cellvec(op,uh)
  assemble(op.assem, cellvec)
end

function apply!(b::AbstractVector,op::NonLinearFEOperator,uh::FEFunction)
  cellvec = _cellvec(op,uh)
  assemble!(b,op.assem, cellvec)
end

function jacobian(op::NonLinearFEOperator,uh::FEFunction)
  cellmat = _cellmat(op,uh)
  assemble(op.assem, cellmat)
end

function jacobian!(mat::AbstractMatrix,op::NonLinearFEOperator,uh::FEFunction)
  cellmat = _cellmat(op,uh)
  assemble!(mat,op.assem, cellmat)
end

function _cellvec(op,uh)
  v = CellBasis(op.testfesp)
  integrate(op.res(uh,v), op.trian, op.quad)
end

function _cellmat(op,uh)
  v = CellBasis(op.testfesp)
  du = CellBasis(op.trialfesp)
  integrate(op.jac(uh,v,du), op.trian, op.quad)
end

struct NonLinearFESolver <: FESolver
  ls::LinearSolver
  tol::Float64
  max_nliters::Int
end

LinearSolver(s::NonLinearFESolver) = s.ls

function solve!(
  uh::FEFunction,A::AbstractMatrix,b::AbstractVector,dx::AbstractVector,ns::NumericalSetup,
  nls::NonLinearFESolver,op::FEOperator)

  # Get raw solution
  x = free_dofs(uh)

  # Check for convergence on the initial residual
  isconv, conv0 = _check_convergence(nls, b)
  if isconv; return; end

  # Newton-like iterations
  for nliter in 1:nls.max_nliters

    # Solve linearized problem
    broadcast!(*,b,b,-1)
    solve!(dx,ns,A,b)
    broadcast!(+,x,x,dx)

    # Check convergence for the current residual
    apply!(b, op, uh)
    isconv = _check_convergence(nls, b, conv0)
    if isconv; return; end

    if nliter == nls.max_nliters
      @unreachable
    end

    # Assemble jacobian (fast in-place version)
    # and prepare solver
    jacobian!(A, op, uh)
    numerical_setup!(ns,A)

  end

end

function _check_convergence(nls::NonLinearFESolver,b)
  m0 = _inf_norm(b)
  (false, m0)
end

function _check_convergence(nls::NonLinearFESolver,b,m0)
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

end # module FEOperators
