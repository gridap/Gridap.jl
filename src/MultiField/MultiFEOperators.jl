include("MultiFEFunctions.jl")

module MultiFEOperators

using Gridap
using Gridap.Helpers
using Gridap.Geometry
using Gridap.CellMaps
using Gridap.CellQuadratures
using Gridap.CellIntegration
using Gridap.LinearSolvers
using Gridap.FESpaces
using Gridap.FEOperators: LinearFEOperator
using Gridap.FEOperators: LinearFESolver
using Gridap.MultiCellArrays
using Gridap.Assemblers
using ..MultiFESpaces
using ..MultiAssemblers
using ..MultiFEFunctions
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

function LinearFEOperator(
  biform::Function,
  liform::Function,
  testfesp::Vector{<:FESpaceWithDirichletData},
  trialfesp::Vector{<:FESpaceWithDirichletData},
  assem::MultiAssembler,
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where Z
  V = MultiFESpace(testfesp)
  U = MultiFESpace(trialfesp)
  LinearMultiFEOperator(biform,liform,V,U,assem,trian,quad)
end

function LinearMultiFEOperator(
  biform::Function,
  liform::Function,
  testfesp::MultiFESpace,
  trialfesp::MultiFESpace,
  assem::MultiAssembler,
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where Z

  # This will not be a CellBasis in the future
  v = [ CellBasis(V) for V in testfesp  ]
  u = [ CellBasis(U) for U in trialfesp ]

  preblocks, fieldids = biform(v,u)

  blocks = [ integrate(bi,trian,quad) for bi in preblocks ]

  cellmat = MultiCellMatrix(blocks,fieldids)

  # The way we modify the rhs can be improved
  uhd = zero(trialfesp)
  preblocks, prefieldids = biform(v,uhd)

  blocks1 = [ integrate(-bi,trian,quad) for bi in preblocks ]
  fieldids1 = [ (fi[1],) for fi in prefieldids ]

  preblocks, fieldids2 = liform(v)

  blocks2 = [ integrate(bi,trian,quad) for bi in preblocks ]

  blocks = vcat(blocks1,blocks2)
  fieldids = vcat(fieldids1,fieldids2)

  cellvec = MultiCellVector(blocks,fieldids)

  mat = assemble(assem,cellmat)
  vec = assemble(assem,cellvec)

  LinearMultiFEOperator(mat,vec,trialfesp,testfesp)

end

TrialFESpace(op::LinearMultiFEOperator) = op.trialfesp

TestFESpace(op::LinearMultiFEOperator) = op.testfesp

function apply(o::LinearMultiFEOperator,uh::MultiFEFunction)
  vals = free_dofs(uh)
  o.mat * vals - o.vec
end

function apply!(b::Vector,o::LinearMultiFEOperator,uh::MultiFEFunction)
  vals = free_dofs(uh)
  mul!(b,o.mat,vals)
  broadcast!(-,b,b,o.vec)
end

jacobian(o::LinearMultiFEOperator,::MultiFEFunction) = o.mat

function jacobian!(
  mat::AbstractMatrix, o::LinearMultiFEOperator, ::MultiFEFunction)
  @assert mat === o.mat
end

struct LinearMultiFESolver <: MultiFESolver
  ls::LinearSolver
end

function solve!(::MultiFEFunction,::LinearMultiFESolver,::MultiFEOperator)
  @unreachable
end

function solve!(::MultiFEFunction,::LinearMultiFESolver,::MultiFEOperator,::Any)
  @unreachable
end

function solve!(uh::MultiFEFunction,s::LinearMultiFESolver,o::LinearMultiFEOperator)
  x = free_dofs(uh)
  A = o.mat
  b = o.vec
  ss = symbolic_setup(s.ls,A)
  ns = numerical_setup(ss,A)
  solve!(x,ns,A,b)
  ns
end

function solve!(uh::MultiFEFunction,s::LinearMultiFESolver,o::LinearMultiFEOperator,ns::NumericalSetup)
  x = free_dofs(uh)
  A = o.mat
  b = o.vec
  solve!(x,ns,A,b)
end

end # module MultiFEOperators
