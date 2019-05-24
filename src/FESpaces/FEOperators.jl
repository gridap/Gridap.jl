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

"""
Solve with a given initial guess that will be overwritten with the solution
"""
function solve!(::FEFunction,::FESolver,::FEOperator)
  @abstractmethod
end

"""
Solve that allocates and returns the solution
"""
function solve(nls::FESolver,op::FEOperator)
  U = TrialFESpace(op)
  uh = zero(U)
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

function solve!(::FEFunction,::LinearFESolver,::FEOperator)
  @unreachable
end

function solve!(uh::FEFunction,s::LinearFESolver,o::LinearFEOperator)
  ss = symbolic_setup(s.ls,o.mat)
  ns = numerical_setup(ss,o.mat)
  x = free_dofs(uh)
  solve!(x,ns,o.mat,o.vec)
  diri_vals = diri_dofs(o.trialfesp)
  FEFunction(o.trialfesp,x,diri_vals)
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

function solve!(uh::FEFunction,nls::NonLinearFESolver,op::FEOperator)

  # Prepare raw initial guess and correction
  x = free_dofs(uh)
  dx = similar(x)

  # Check convergence for the initial residual
  b = apply(op, uh)
  isconv, conv0 = _check_convergence(nls, b)
  if isconv; return; end

  # Assemble Jacobian and prepare solver
  A = jacobian(op, uh)
  ss = symbolic_setup(nls.ls, A)
  ns = numerical_setup(ss,A)

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
    numerical_setup!(ns,ss,A)

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
