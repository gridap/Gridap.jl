module FEOperators

using Gridap.Helpers
using Gridap.Geometry
using Gridap.CellMaps
using Gridap.CellQuadratures
using Gridap.CellIntegration

using ..FESpaces
using ..Assemblers

export FEOperator
export FESolver
export LinearFEOperator
export LinearFESolver
export solve

abstract type FEOperator end

function apply(::FEOperator,::FEFunction)::AbstractVector
  @abstractmethod
end

function jacobian(::FEOperator,::FEFunction)::AbstractMatrix
  @abstractmethod
end

abstract type FESolver end

function solve(::FESolver,::FEOperator)::FEFunction
  @abstractmethod
end

"""
Struct representing a linear FE Operator
"""
struct LinearFEOperator{M,V} <:FEOperator
  mat::M
  vec::V
  trialfesp::FESpaceWithDirichletData
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

  cellmat = integrate(biform(v,u),trian,quad)

  # The way we modify the rhs can be improved
  fun(x) = zero(T)
  uhd = interpolate(trialfesp,fun)

  cellvec = integrate( liform(v)-biform(v,uhd), trian, quad)

  mat = assemble(assem,cellmat)

  vec = assemble(assem,cellvec)

  LinearFEOperator(mat,vec,trialfesp)

end

function apply(o::LinearFEOperator,uh::FEFunction)::AbstractVector
  vals = free_dofs(uh)
  o.mat * vals - o.vec
end

jacobian(o::LinearFEOperator,::FEFunction) = o.mat

struct LinearFESolver <: FESolver end

function solve(::LinearFESolver,::FEOperator)
  @unreachable
end

function solve(s::LinearFESolver,o::LinearFEOperator)
  x = o.mat \ o.vec
  diri_vals = diri_dofs(o.trialfesp)
  FEFunction(o.trialfesp,x,diri_vals)
end

end # module FEOperators
