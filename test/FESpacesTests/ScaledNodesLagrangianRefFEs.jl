# module ModifiedNodesLagRefFEs
#
# using Gridap


function ScaledNodesLagrangianRefFE(::Type{T}, p::Polytope{D}, orders::Vector{Int}, scaling) where {D,T}
  @assert length(orders) == D
  @assert D > 0
  if !(all(extrusion(p).array.== HEX_AXIS) || all(orders.==orders[1]))
    error("One can consider anisotropic orders on n-cubes only")
  end
  if (all(orders.==1))
    nodes, nfacenodes = Gridap.RefFEs._linear_lagrangian_nodes_polytope(p)
  else
    nodes, nfacenodes = Gridap.RefFEs._high_order_lagrangian_nodes_polytope(p,orders)
  end
  nodes = nodes*scaling
  dofsb = Gridap.RefFEs.LagrangianDOFBasis(T,nodes)
  # @santiagobadia : This part is missing for anisotropic orders !
  prebasis = Gridap.RefFEs._monomial_basis(p,T,orders[1])
  aux = zeros(Float64,numlocaldofs(dofsb),numlocaldofs(dofsb))
  @assert numlocaldofs(dofsb) == length(prebasis)
  evaluate!(dofsb,prebasis,aux)
  changeofbasis=inv(aux)
  basis = change_basis(prebasis, changeofbasis)
  nfacedofs = Gridap.RefFEs._genereate_nface_to_dofs(nfacenodes,dofsb.node_and_comp_to_dof)
  return LagrangianRefFE{D,T}(p, dofsb, basis, nfacedofs, nfacenodes)
end

function ScaledNodesLagrangianRefFE(::Type{T}, p::Polytope{D}, order::Int, scaling) where {D,T}
  _order = order*ones(Int,D)
  return ScaledNodesLagrangianRefFE(T,p,_order,scaling)
end

# F = Gridap.RefFEs
#
# T = Float64
# D = 1
# model = CartesianDiscreteModel(partition=(1,))
# diri = [1]
# p = Polytope(1,)
# order = 2
# reffe1 = LagrangianRefFE(Float64,p,order)
# orders = order*ones(Int,D)
# reffe2 = ScaledLagrangianRefFE(Float64,p,orders,0.5)
#
#
# end # module
