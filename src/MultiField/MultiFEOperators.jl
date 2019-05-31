include("MultiFEFunctions.jl")

module MultiFEOperators

using Gridap
using Gridap.Helpers
using Gridap.FESpaces
using ..MultiFESpaces
using ..MultiAssemblers
using ..MultiFEFunctions

abstract type MultiFEOperator end

function apply(::MultiFEOperator,::MultiFEFunction)::AbstractVector
  @abstractmethod
end

function jacobian(::MultiFEOperator,::MultiFEFunction)::AbstractMatrix
  @abstractmethod
end

function apply!(::AbstractVector,::MultiFEOperator,::MultiFEFunction)
  @abstractmethod
end

function jacobian!(::AbstractMatrix,::MultiFEOperator,::MultiFEFunction)
  @abstractmethod
end

function TrialFESpace(::MultiFEOperator)::FESpaceWithDirichletData
  @abstractmethod
end

function TestFESpace(::MultiFEOperator)::FESpaceWithDirichletData
  @abstractmethod
end

function residual(op::MultiFEOperator,uh::MultiFEFunction)
  apply(op,uh)
end

function residual!(b::AbstractVector,op::MultiFEOperator,uh::MultiFEFunction)
  apply!(b,op,uh)
end

abstract type MultiFESolver end

function solve!(uh::MultiFEFunction,::MultiFESolver,::MultiFEOperator)::Any
  @abstractmethod
end

function solve!(uh::MultiFEFunction,::MultiFESolver,::MultiFEOperator,::Any)
  @abstractmethod
end

function solve(nls::MultiFESolver,op::MultiFEOperator)
  U = TrialFESpace(op)
  uh = zero(U)
  solve!(uh,nls,op)
  uh
end

struct LinearMultiFEOperator{M,V} <:MultiFEOperator
  mat::M
  vec::V
  trialfesp::MultiFESpace
  testfesp::MultiFESpace
end

function LinearMultiFEOperator(
  biform::Function,
  liform::Function,
  testfesp::MultiFESpace,
  trialfesp::MultiFESpace,
  assem::Assembler{M,V},
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where {M,V,Z}

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

end # module MultiFEOperators
