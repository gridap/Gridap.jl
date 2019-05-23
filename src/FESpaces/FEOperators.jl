module FEOperators

using Gridap.Helpers
using Gridap.Geometry
using Gridap.CellValues
using Gridap.CellMaps
using Gridap.CellQuadratures
using Gridap.CellIntegration

using Gridap.FESpaces
using Gridap.Assemblers

export FEOperator
export FESolver
export LinearFEOperator
export LinearFESolver
export NonLinearFEOperator
export NonLinearFESolver
export solve
export jacobian
import Gridap: apply

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

  # The way we modify the rhs can be improved
  E = eltype(T)
  free_values = zeros(E,num_free_dofs(trialfesp))
  diri_values = diri_dofs(trialfesp)
  uhd = FEFunction(trialfesp,free_values,diri_values)

  cellmat = integrate(biform(v,u),trian,quad)
  cellvec = integrate( liform(v)-biform(v,uhd), trian, quad)

  mat = assemble(assem,cellmat)
  vec = assemble(assem,cellvec)

  LinearFEOperator(mat,vec,trialfesp)

end

function apply(o::LinearFEOperator,uh::FEFunction)
  vals = free_dofs(uh)
  o.mat * vals - o.vec
end

jacobian(o::LinearFEOperator,::FEFunction) = o.mat

"""
The solver that solves a LinearFEOperator
"""
struct LinearFESolver <: FESolver end

function solve(::LinearFESolver,::FEOperator)
  @unreachable
end

function solve(s::LinearFESolver,o::LinearFEOperator)
  x = o.mat \ o.vec
  diri_vals = diri_dofs(o.trialfesp)
  FEFunction(o.trialfesp,x,diri_vals)
end

"""
Struct representing a non-linear FE Operator
"""
struct NonLinearFEOperator{M,V,E} <:FEOperator
  uh::FEFunction
  cellvec::CellVector{E}
  cellmat::CellMatrix{E}
  assem::Assembler{M,V}
end

function NonLinearFEOperator(
  res::Function,
  jac::Function,
  testfesp::FESpace{D,Z,T},
  trialfesp::FESpaceWithDirichletData,
  assem::Assembler{M,V},
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where {M,V,D,Z,T}

  # Construction of the initial guess
  E = eltype(T)
  free_values = zeros(E,num_free_dofs(trialfesp))
  diri_values = diri_dofs(trialfesp)
  uh = FEFunction(trialfesp,free_values,diri_values)

  # This will not be a CellBasis in the future
  v = CellBasis(testfesp)
  du = CellBasis(trialfesp)

  cellmat = integrate(jac(uh,v,du), trian, quad)
  cellvec = integrate(res(uh,v), trian, quad)

  NonLinearFEOperator(uh,cellvec,cellmat,assem)

end

function apply(op::NonLinearFEOperator,uh::FEFunction)
  _update_free_values(op,uh)
  assemble(op.assem, op.cellvec)
end

function jacobian(op::NonLinearFEOperator,uh::FEFunction)
  _update_free_values(op,uh)
  assemble(op.assem, op.cellmat)
end

function _update_free_values(op,uh)
  op_free_values = free_dofs(op.uh)
  free_values = free_dofs(uh)
  if !(op_free_values === free_values)
    op_free_values .= free_values
  end
end

"""
The solver that solves a NonLinearFEOperator
"""
struct NonLinearFESolver <: FESolver
  tol::Float64
  maxiters::Int
end

function solve(::NonLinearFESolver,::FEOperator)
  @unreachable
end

# For the moment it can only solve a NonLinearFEOperator
function solve(s::NonLinearFESolver,o::NonLinearFEOperator)

  # Get the solution vector
  uh = o.uh
  x = free_dofs(uh)

  # Prepare initial guess
  T = value_type(uh)
  E = eltype(T)
  x .= zero(E)

  # For, efficiency we can introduce in place variants apply! and jacobian!
  
  b = apply(o,uh)
  m0 = maximum(abs.(b))

  max_nliters = s.maxiters
  tol = s.tol

  for nliter in 1:max_nliters

    A = jacobian(o,uh)
    b *= -1
    dx = A \ b
    broadcast!(+,x,x,dx)

    b .= apply(o,uh)

    m = maximum(abs.(b))
    if  m < tol * m0
      break
    end

    if nliter == max_nliters
      @unreachable
    end

  end

  uh

end

end # module FEOperators
