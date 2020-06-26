function CDLagrangianRefFE(::Type{T},p::Polytope{D},order::Int) where {T,D}
  orders = tfill(order,Val{D}())
  CDLagrangianRefFE(T,p,orders)
end

function CDLagrangianRefFE(::Type{T},p::Polytope{D},orders,cd_dim) where {T,D}
  lrfe = LagrangianRefFE(T,p,orders)
  grfe = lrfe.data

  # cd_dim could be a vector that tells you whether the dim is
  # CG or DG and use it to build face_xxx


  # Now change face_own_dofs
  # face_own_dofs_permutations can be an identity map
  # face_dofs are the nodes in the clausure of a face, e.g.,
  # the one of the face and the faces in their boundary,
  # e.g., for a face, the dofs owned by the face and its
  # vertices and edges on the face boundary
  # probably take a look at the function
  # _generate_face_dofs(ndofs,face_to_own_dofs,polytope,reffaces)
  GenericRefFE{D}(
      grfe.ndofs,
      grfe.polytope,
      grfe.prebasis,
      grfe.dofs,
      face_own_dofs,
      face_own_dofs_permutations,
      face_dofs,
      grfe.shapefuns)
end


#
# function TensorProductReferenceFE(reffe1::ReferenceFE, reffe2::ReferenceFE)
#   ...
#   GenericRefFE
# end
#
# function TensorProductReferenceFE(reffe1::NodalReferenceFE, reffe2::NodalReferenceFE)
#   ...
#   GenericNodalCartesianReferenceFE
# end
#
# function TensorProductReferenceFE(reffe1::NodalReferenceFE, reffe2::NodalReferenceFE)
#   @abstractmethod
# end
#
# # struct DiscRefFE{D} <: NodalReferenceFE{D}
# #   p_reffe::LagrangianRefFE{D}
# #   polytope::Polytope{D}
# # end
# #
# # function PDiscRefFE(::Type{T},p::Polytope,order::Integer) where T
# #   D = num_cell_dims(p)
# #   extrusion = tfill(TET_AXIS,Val{D}())
# #   simplex = ExtrusionPolytope(extrusion)
# #   p_reffe = LagrangianRefFE(T,simplex,order)
# #   new{D}(p_reffe,p)
# # end
#
# function QDiscRefFE(::Type{T},p::Polytope,order::Integer) where T
#   D = num_cell_dims(p)
#   extrusion = tfill(HEX_AXIS,Val{D}())
#   hex = ExtrusionPolytope(extrusion)
#   p_reffe = LagrangianRefFE(T,hex,order)
#   PDiscRefFE{D}(p_reffe,p)
# end
#
# TensorProductReferenceFE(reffe1::PDiscRefFE, reffe2::LagrangianRefFE) = TensorProductReferenceFE(reffe2,reffe1)
#
# function TensorProductReferenceFE(reffe1::LagrangianRefFE, reffe2::PDiscRefFE)
#
#   T, p, orders =  _tensor_product_reference_fe(reffe1,reffe2)
#
#   is_hex(p1), is_hex(p2)
#
#   prebasis = compute_monomial_basis(T,p,orders)
#   nodes, face_own_nodes = compute_nodes(p,orders)
#   # face_own_nodes = tensor_product_face_own_nodes(fon1,fon2)
#   dofs = LagrangianDofBasis(T,nodes)
#   interior_nodes = dofs.nodes[face_own_nodes[end]]
#   own_nodes_permutations = [ Int[] ]
#   face_own_nodes = [ Int[] ]
#   # own_nodes_permutations = compute_own_nodes_permutations(p, interior_nodes)
#   # reffaces = compute_tensor_product_lagrangian_reffaces(T,p,orders)
#   LagrangianRefFE(p,prebasis,dofs,face_own_nodes,own_nodes_permutations,reffaces)
# end
#
# function TensorProductReferenceFE(reffe1::LagrangianRefFE, reffe2::LagrangianRefFE)
#   T, p, order =  _tensor_product_reference_fe(reffe1,reffe2)
#   LagrangianRefFE(T1,p,order)
# end
#
# function TensorProductReferenceFE(reffe1::LagrangianRefFE, reffe2::LagrangianRefFE)
#   T, p, order =  _tensor_product_reference_fe(reffe1,reffe2)
#   LagrangianRefFE(T1,p,order)
# end
#
# function TensorProductReferenceFE(reffe1::PDiscRefFE, reffe2::PDiscRefFE)
#   T, p, order =  _tensor_product_reference_fe(reffe1,reffe2)
#   PDiscRefFE(T1,p,order)
# end
#
# function _tensor_product_reference_fe(reffe1,reffe1)
#
#   p1 = get_polytope(reffe1)
#   p2 = get_polytope(reffe2)
#
#   e1 = get_extrusion(p1)
#   e2 = get_extrusion(p2)
#
#   e = (e1...,e2...)
#
#   p = Polytope(e)
#
#   order1 = get_order(reffe1)
#   order2 = get_order(reffe2)
#
#   order = (order1..., order2...)
#
#   pb1 = get_prebasis(reffe1)
#   pb2 = get_prebasis(reffe2)
#
#   T1 = get_value_type(pb1)
#   T2 = get_value_type(pb2)
#
#   # Check whether it can be more general
#   if (T1 != T2)
#     error("Types of ref FEs for tensor product must be the same")
#   end
#
#   return T1, p, order
# end
