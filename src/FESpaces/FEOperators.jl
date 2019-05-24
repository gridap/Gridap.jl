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
import Gridap: solve
export jacobian
import Gridap: apply

abstract type FEOperator end

# @santiagobadia : Public method
function apply(::FEOperator,::FEFunction)::AbstractVector
  @abstractmethod
end

# @santiagobadia : Not sure it has to be public
function jacobian(::FEOperator,::FEFunction)::AbstractMatrix
  @abstractmethod
end

abstract type FESolver end

# @santiagobadia : Public method
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
  # @santiagobadia : CellBases + block-data for multifield problems
  v = CellBasis(testfesp)
  u = CellBasis(trialfesp)

  # The way we modify the rhs can be improved
  E = eltype(T)
  free_values = zeros(E,num_free_dofs(trialfesp))
  diri_values = diri_dofs(trialfesp)
  uhd = FEFunction(trialfesp,free_values,diri_values)

  cellmat = integrate(biform(v,u),trian,quad)
  cellvec = integrate( liform(v)-biform(v,uhd), trian, quad)
  # @santiagobadia : Where is the problem with the RHS computation?
  # It is nice and I don't see any efficiency issue

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
  @unreachable # @santiagobadia : Why unreachable? Not implemented?
end

function solve(s::LinearFESolver,o::LinearFEOperator)
  x = o.mat \ o.vec
  diri_vals = diri_dofs(o.trialfesp)
  FEFunction(o.trialfesp,x,diri_vals)
end

"""
Struct representing a nonlinear FE Operator
"""
struct NonLinearFEOperator{M,V,E} <:FEOperator
  uh::FEFunction
  # @santiagobadia : The cellvec and cellmat are needed just as scratch data
  # for efficiency
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
  # @santiagobadia : Why are we putting to zero x, probably we can use uh as
  # initial guess

  # For, efficiency we can introduce in place variants apply! and jacobian!

  # @santiagobadia:
  # I would say that all the lines below should go the NonlinearFESolver
  # since different implementations are needed (relaxation, Anderson acc,
  # nonlinear GMRES, etc), different norms, etc. In any case, we can probably
  # define an API for for nonlinear solvers such that we can work a
  # generic nonlinear loop...
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
      @unreachable # @santiagobadia : @notconverged better?
    end

  end

  uh

end

end # module FEOperators
