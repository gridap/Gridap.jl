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

abstract type FESolver end

function solve!(uh::MultiFEFunction,::FESolver,::MultiFEOperator)::Any
  @abstractmethod
end

function solve!(uh::MultiFEFunction,::FESolver,::MultiFEOperator,::Any)
  @abstractmethod
end

function solve(nls::FESolver,op::MultiFEOperator)
  U = TrialFESpace(op)
  uh = zero(U)
  solve!(uh,nls,op)
  uh
end

end # module MultiFEOperators
