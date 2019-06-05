include("MultiFEBases.jl")

module MultiFEOperators

using Gridap
using Gridap.Geometry
using Gridap.CellQuadratures
using Gridap.FESpaces
using Gridap.FEOperators: LinearFEOperator
using Gridap.FEOperators: NonLinearFEOperator
using ..MultiFESpaces
using ..MultiAssemblers
using ..MultiFEFunctions
using ..MultiFEBases

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
  LinearFEOperator(biform,liform,V,U,assem,trian,quad)
end

function NonLinearFEOperator(
  res::Function,
  jac::Function,
  testfesp::Vector{<:FESpaceWithDirichletData},
  trialfesp::Vector{<:FESpaceWithDirichletData},
  assem::MultiAssembler,
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where Z
  V = MultiFESpace(testfesp)
  U = MultiFESpace(trialfesp)
  NonLinearFEOperator(res,jac,V,U,assem,trian,quad)
end

end # module MultiFEOperators
