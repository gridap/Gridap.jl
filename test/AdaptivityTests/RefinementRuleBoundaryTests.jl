module RefinementRuleBoundaryTests

using Test
using Gridap
using Gridap.Helpers
using Gridap.Adaptivity
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays

function _get_terms(poly::Polytope,orders)
  _nodes, facenodes = Gridap.ReferenceFEs._compute_nodes(poly,orders)
  terms = Gridap.ReferenceFEs._coords_to_terms(_nodes,orders)
  return terms
end

function _get_face_orders(p::Polytope{Dc},D::Int,orders::Tuple) where Dc
  @check length(orders) == Dc
  @check 1 <= D < Dc
  @check is_n_cube(p)

  if D == 1 # Edges (2D, 3D)
    tangents = get_edge_tangent(p)
    face_orders = map(tangents) do t
      axis = findfirst(i -> abs(t[i]) > 0.5 ,1:Dc)
      return [orders[axis]]
    end
  elseif D == Dc-1 # Faces (3D)
    normals = get_facet_normal(p)
    face_orders = map(normals) do n
      mask = map(i -> abs(n[i]) < 1.e-3,1:Dc)
      return [orders[mask]...]
    end
  else
    @notimplemented
  end

  return face_orders
end

D = 1
rr = Gridap.Adaptivity.RedRefinementRule(QUAD)
poly  = get_polytope(rr)
coarse_orders = (4,4)
coarse_reffe  = ReferenceFE(poly,lagrangian,Float64,coarse_orders)
coarse_face_nodes = coarse_reffe.face_nodes
coarse_face_polys = CompressedArray(Gridap.ReferenceFEs._compute_reffaces_and_face_types(poly,Val(D))...)
c_edge_to_coarse_dof = coarse_face_nodes[get_dimranges(poly)[D+1]]

cell_polys  = Gridap.Adaptivity.get_cell_polytopes(rr)
fine_orders = coarse_orders .÷ 2
fine_reffe  = lazy_map(p->ReferenceFE(p,lagrangian,Float64,fine_orders),cell_polys)

model = Gridap.Adaptivity.get_ref_grid(rr)

fine_topo = get_grid_topology(model)
fine_boundary_faces = findall(get_isboundary_face(fine_topo,D))
fine_face_grid = Grid(ReferenceFE{D},model)
fine_face_polys = CompressedArray(map(get_polytope,get_reffes(fine_face_grid)),get_cell_type(fine_face_grid))
fine_boundary_polys = lazy_map(Reindex(fine_face_polys),fine_boundary_faces)

d_to_face_to_child_faces = Gridap.Adaptivity.get_d_to_face_to_child_faces(rr)
face_to_child_faces = d_to_face_to_child_faces[D+1]

coarse_face_orders = _get_face_orders(poly,D,coarse_orders)
fine_face_orders = _get_face_orders(poly,D,fine_orders)

num_coarse_faces = num_faces(coarse_reffe,D)
coarse_dofs_above_fine_dofs = Vector{Vector{Vector{Int32}}}(undef,num_coarse_faces)
for cF in 1:num_coarse_faces
  coarse_face_poly = coarse_face_polys[cF]
  coarse_terms = _get_terms(coarse_face_poly,coarse_face_orders[cF])
  coarse_dofs  = zeros(Int32,Tuple(maximum(coarse_terms)))
  coarse_dofs[coarse_terms] .= c_edge_to_coarse_dof[cF]

  child_faces = face_to_child_faces[cF]
  fine_dofs = Vector{Vector{Int32}}(undef,length(child_faces))
  for (i,fF) in enumerate(child_faces)
    fine_face_poly = fine_face_polys[fF]
    fine_terms = _get_terms(fine_face_poly,fine_face_orders[cF])

    local_dof_range = map(o->(i-1)*o+1:i*o+1,fine_face_orders[cF])
    local_coarse_dofs = view(coarse_dofs,local_dof_range...)
    fine_dofs[i] = map(Reindex(local_coarse_dofs),fine_terms)
  end
  coarse_dofs_above_fine_dofs[cF] = fine_dofs
end


# Test the final function

res = Adaptivity.get_face_subface_ldof_to_cell_ldof(rr,(2,2),1)

end