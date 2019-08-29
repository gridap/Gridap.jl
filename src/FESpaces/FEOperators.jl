module FEOperators

using Gridap
using Gridap.Helpers

using LinearAlgebra

import Gridap: apply
export apply!
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
import Gridap: residual
import Gridap: residual!
import Gridap.NonLinearSolvers: create_in_domain
import Gridap.FESpaces: TrialFESpace
import Gridap.FESpaces: TestFESpace
import Gridap.LinearSolvers: LinearSolver

# We use duck typing in this module for the following types:

"""
Arguments annotated with this type have to implement the following queries:

free_dofs(::FEFunctionLike)

"""
const FEFunctionLike = Any

"""
Arguments annotated with this type have to implement the following queries:

FEFunction(::FESpaceLike, free_vals)

FEBasis(::FESpaceLike)

zero(::FESpaceLike)

"""
const FESpaceLike = Any

"""
Arguments annotated with this type have to implement the following queries:

assemble(::AssemblerLike, cellmat)

assemble(::AssemblerLike, cellvec)

assemble!(A,::AssemblerLike, cellmat)

assemble!(b,::AssemblerLike, cellvec)

"""
const AssemblerLike = Any

abstract type FEOperator end

function apply(::FEOperator,::FEFunctionLike)::AbstractVector
  @abstractmethod
end

function jacobian(::FEOperator,::FEFunctionLike)::AbstractMatrix
  @abstractmethod
end

function apply!(::AbstractVector,::FEOperator,::FEFunctionLike)
  @abstractmethod
end

function jacobian!(::AbstractMatrix,::FEOperator,::FEFunctionLike)
  @abstractmethod
end

function TrialFESpace(::FEOperator)::FESpaceLike
  @abstractmethod
end

function TestFESpace(::FEOperator)::FESpaceLike
  @abstractmethod
end

function residual(op::FEOperator,uh::FEFunctionLike)
  apply(op,uh)
end

function residual!(b::AbstractVector,op::FEOperator,uh::FEFunctionLike)
  apply!(b,op,uh)
end

abstract type FESolver end

function solve!(uh::FEFunctionLike,::FESolver,::FEOperator)::Any
  @abstractmethod
end

function solve!(uh::FEFunctionLike,::FESolver,::FEOperator,::Any)
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
  zh = zero(U)
  free_dofs(zh)
end

"""
Struct representing a linear FE Operator
"""
struct LinearFEOperator{M,V} <:FEOperator
  mat::M
  vec::V
  trialfesp::FESpaceLike
  testfesp::FESpaceLike
end

function LinearFEOperator(
  testfesp::FESpaceLike,
  trialfesp::FESpaceLike,
  assem::AssemblerLike,
  terms::Vararg{<:AffineFETerm})

  v = FEBasis(testfesp)
  u = FEBasis(trialfesp)

  uhd = zero(trialfesp)

  cms = setup_cell_matrix(v,u,terms...)
  cvs = setup_cell_vector(v,uhd,terms...)

  mat = assemble(assem,cms...)
  vec = assemble(assem,cvs...)

  LinearFEOperator(mat,vec,trialfesp,testfesp)

end

function LinearFEOperator(
  testfesp::FESpaceLike,
  trialfesp::FESpaceLike,
  terms::Vararg{<:AffineFETerm})

  assem = SparseMatrixAssembler(testfesp,trialfesp)
  LinearFEOperator(testfesp,trialfesp,assem,terms...)
end

function LinearFEOperator(
  biform::Function,
  liform::Function,
  testfesp::FESpaceLike,
  trialfesp::FESpaceLike,
  assem::AssemblerLike,
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where Z

  term = AffineFETerm(biform,liform,trian,quad)

  LinearFEOperator(testfesp,trialfesp,assem,term)

end

TrialFESpace(op::LinearFEOperator) = op.trialfesp

TestFESpace(op::LinearFEOperator) = op.testfesp

function apply(o::LinearFEOperator,uh::FEFunctionLike)
  vals = free_dofs(uh)
  o.mat * vals - o.vec
end

function apply!(b::Vector,o::LinearFEOperator,uh::FEFunctionLike)
  vals = free_dofs(uh)
  mul!(b,o.mat,vals)
  broadcast!(-,b,b,o.vec)
end

jacobian(o::LinearFEOperator,::FEFunctionLike) = o.mat

function jacobian!(mat::AbstractMatrix, o::LinearFEOperator, ::FEFunctionLike)
  @assert mat === o.mat
end

"""
The solver that solves a LinearFEOperator
"""
struct LinearFESolver <: FESolver
  ls::LinearSolver
end

function solve!(::FEFunctionLike,::LinearFESolver,::FEOperator)
  @unreachable
end

function solve!(::FEFunctionLike,::LinearFESolver,::FEOperator,::Any)
  @unreachable
end

function solve!(uh::FEFunctionLike,s::LinearFESolver,o::LinearFEOperator)
  x = free_dofs(uh)
  A = o.mat
  b = o.vec
  ss = symbolic_setup(s.ls,A)
  ns = numerical_setup(ss,A)
  solve!(x,ns,A,b)
  ns
end

function solve!(uh::FEFunctionLike,s::LinearFESolver,o::LinearFEOperator,ns::NumericalSetup)
  x = free_dofs(uh)
  A = o.mat
  b = o.vec
  solve!(x,ns,A,b)
end

function solve(op::LinearFEOperator)
  ls = LUSolver()
  solver = LinearFESolver(ls)
  solve(solver,op)
end

"""
Struct representing a nonlinear FE Operator
"""
struct NonLinearFEOperator <:FEOperator
  testfesp::FESpaceLike
  trialfesp::FESpaceLike
  assem::AssemblerLike
  terms::NTuple{N,<:FETerm} where N
end

function NonLinearFEOperator(
  testfesp::FESpaceLike,
  trialfesp::FESpaceLike,
  assem::AssemblerLike,
  terms::Vararg{<:FETerm})
  NonLinearFEOperator(testfesp,trialfesp,assem,terms)
end

function NonLinearFEOperator(
  testfesp::FESpaceLike,
  trialfesp::FESpaceLike,
  terms::Vararg{<:FETerm})
  assem = SparseMatrixAssembler(testfesp,trialfesp)
  NonLinearFEOperator(testfesp,trialfesp,assem,terms)
end

function NonLinearFEOperator(
  res::Function,
  jac::Function,
  testfesp::FESpaceLike,
  trialfesp::FESpaceLike,
  assem::AssemblerLike,
  trian::Triangulation,
  quad::CellQuadrature)

  term = NonLinearFETerm(res,jac,trian,quad)

  NonLinearFEOperator(testfesp,trialfesp,assem,term)

end

TrialFESpace(op::NonLinearFEOperator) = op.trialfesp

TestFESpace(op::NonLinearFEOperator) = op.testfesp

function apply(op::NonLinearFEOperator,uh::FEFunctionLike)
  cellvec = _cellvec(op,uh)
  assemble(op.assem, cellvec...)
end

function apply!(b::AbstractVector,op::NonLinearFEOperator,uh::FEFunctionLike)
  cellvec = _cellvec(op,uh)
  assemble!(b,op.assem, cellvec...)
end

function jacobian(op::NonLinearFEOperator,uh::FEFunctionLike)
  cellmat = _cellmat(op,uh)
  assemble(op.assem, cellmat...)
end

function jacobian!(mat::AbstractMatrix,op::NonLinearFEOperator,uh::FEFunctionLike)
  cellmat = _cellmat(op,uh)
  assemble!(mat,op.assem, cellmat...)
end

function _cellvec(op,uh)
  v = FEBasis(op.testfesp)
  setup_cell_residual(uh,v,op.terms...)
end

function _cellmat(op,uh)
  v = FEBasis(op.testfesp)
  du = FEBasis(op.trialfesp)
  setup_cell_jacobian(uh,v,du,op.terms...)
end

"""
A general NonLinearFESolver
"""
struct NonLinearFESolver <: FESolver
  nls::NonLinearSolver
end

function solve!(uh::FEFunctionLike,nls::NonLinearFESolver,op::FEOperator)
  nlop = NonLinearOpFromFEOp(op)
  x = free_dofs(uh)
  solve!(x,nls.nls,nlop)
end

function solve!(uh::FEFunctionLike,nls::NonLinearFESolver,op::FEOperator,cache::Any)
  nlop = NonLinearOpFromFEOp(op)
  x = free_dofs(uh)
  solve!(x,nls.nls,nlop,cache)
end

end # module FEOperators
