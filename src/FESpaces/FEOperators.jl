module FEOperators

using Gridap.Helpers
using Gridap.Geometry
using Gridap.CellValues
using Gridap.CellMaps
using Gridap.CellQuadratures
using Gridap.CellIntegration

using Gridap.FESpaces
using Gridap.Assemblers
using Gridap.NonLinearSolvers
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
import Gridap: jacobian
import Gridap: jacobian!
import Gridap: apply
import Gridap: apply!
import Gridap: residual
import Gridap: residual!
import Gridap.NonLinearSolvers: create_in_domain
import Gridap.FESpaces: TrialFESpace
import Gridap.FESpaces: TestFESpace
import Gridap.LinearSolvers: LinearSolver

abstract type FEOperator end

# @santiagobadia : Public method
function apply(::FEOperator,::FEFunction)::AbstractVector
  @abstractmethod
end

# @santiagobadia : Not sure it has to be public
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

function residual(op::FEOperator,uh::FEFunction)
  apply(op,uh)
end

function residual!(b::AbstractVector,op::FEOperator,uh::FEFunction)
  apply!(b,op,uh)
end

abstract type FESolver end

function solve!(uh::FEFunction,::FESolver,::FEOperator)::Any
  @abstractmethod
end

function solve!(uh::FEFunction,::FESolver,::FEOperator,::Any)
  @abstractmethod
end

"""
Solve that allocates, and sets initial guess to zero
and returns the solution
"""
function solve(nls::FESolver,op::FEOperator)
  U = TrialFESpace(op)
  uh = zero(U)
  solve!(uh,nls,op)
  uh
end

"""
Struct representing the non-linear algebraic problem
associated with a FEOperator
"""
struct NonLinearOpFromFEOp <: NonLinearOperator
  feop::FEOperator
end

function residual!(b::AbstractVector,op::NonLinearOpFromFEOp,x::AbstractVector)
  U = TrialFESpace(op.feop)
  uh = FEFunction(U,x)
  apply!(b,op.feop,uh)
end

function jacobian!(A::AbstractMatrix,op::NonLinearOpFromFEOp,x::AbstractVector)
  U = TrialFESpace(op.feop)
  uh = FEFunction(U,x)
  jacobian!(A,op.feop,uh)
end

function jacobian(op::NonLinearOpFromFEOp,x::AbstractVector)
  U = TrialFESpace(op.feop)
  uh = FEFunction(U,x)
  jacobian(op.feop,uh)
end

function create_in_domain(op::NonLinearOpFromFEOp)
  U = TrialFESpace(op.feop)
  T = value_type(U)
  E = eltype(T)
  n = length(free_dofs(U))
  x = zeros(E,n)
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

  v = FEBasis(testfesp)
  u = FEBasis(trialfesp)

  uhd = zero(trialfesp)

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

function solve!(::FEFunction,::LinearFESolver,::FEOperator,::Any)
  @unreachable
end

function solve!(uh::FEFunction,s::LinearFESolver,o::LinearFEOperator)
  x = free_dofs(uh)
  A = o.mat
  b = o.vec
  ss = symbolic_setup(s.ls,A)
  ns = numerical_setup(ss,A)
  solve!(x,ns,A,b)
  ns
end

function solve!(uh::FEFunction,s::LinearFESolver,o::LinearFEOperator,ns::NumericalSetup)
  x = free_dofs(uh)
  A = o.mat
  b = o.vec
  solve!(x,ns,A,b)
end

"""
Struct representing a nonlinear FE Operator
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
  v = FEBasis(op.testfesp)
  integrate(op.res(uh,v), op.trian, op.quad)
end

function _cellmat(op,uh)
  v = FEBasis(op.testfesp)
  du = FEBasis(op.trialfesp)
  integrate(op.jac(uh,v,du), op.trian, op.quad)
end


"""
A general NonLinearFESolver
"""
struct NonLinearFESolver <: FESolver
  nls::NonLinearSolver
end

function solve!(uh::FEFunction,nls::NonLinearFESolver,op::FEOperator)
  nlop = NonLinearOpFromFEOp(op)
  x = free_dofs(uh)
  solve!(x,nls.nls,nlop)
end

function solve!(uh::FEFunction,nls::NonLinearFESolver,op::FEOperator,cache::Any)
  nlop = NonLinearOpFromFEOp(op)
  x = free_dofs(uh)
  solve!(x,nls.nls,nlop,cache)
end

end # module FEOperators
