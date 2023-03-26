module FaceLabelingTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.Adaptivity

Dc = 2
model  = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(1,1)))
ctopo = get_grid_topology(model)

fmodel = refine(model)
ftopo  = get_grid_topology(fmodel)

coarse_labeling = get_face_labeling(model)
add_tag_from_tags!(coarse_labeling,"left_boundary",[7])
add_tag_from_tags!(coarse_labeling,"bottom_corners",[1,2])

# d_to_fface_to_cface maps
glue   = Adaptivity.get_adaptivity_glue(fmodel)
rrules = Adaptivity.get_old_cell_refinement_rules(glue)
fine_rrules = Adaptivity.get_new_cell_refinement_rules(glue)

ccell_to_d_to_faces = lazy_map(rr->map(d->Gridap.Geometry.get_faces(get_grid_topology(rr.ref_grid),Dc,d),0:Dc),rrules)
ccell_to_d_to_fface_to_parent_face = lazy_map(Adaptivity.get_d_to_face_to_parent_face,rrules)

ccell_to_fcell = glue.o2n_faces_map
d_to_ccell_to_cface = map(d->Geometry.get_faces(ctopo,Dc,d),0:Dc)
d_to_fcell_to_fface = map(d->Geometry.get_faces(ftopo,Dc,d),0:Dc)

d_to_fface_to_cface = [fill(Int32(0),num_faces(ftopo,d)) for d in 0:Dc]
d_to_fface_to_cface_dim = [fill(Int32(0),num_faces(ftopo,d)) for d in 0:Dc]

for ccell in 1:num_cells(model)
  local_d_to_fface_to_parent_face,
    local_d_to_fface_to_parent_dim = ccell_to_d_to_fface_to_parent_face[ccell]

  for (child,fcell) in enumerate(ccell_to_fcell[ccell])
    for d in 0:Dc
      for (iF,fface) in enumerate(d_to_fcell_to_fface[d+1][fcell])
        fface_coarse_id = ccell_to_d_to_faces[ccell][d+1][child][iF]
        parent = local_d_to_fface_to_parent_face[d+1][fface_coarse_id]
        cface = d_to_ccell_to_cface[d+1][ccell][parent]
        cface_dim = local_d_to_fface_to_parent_dim[d+1][fface_coarse_id]

        d_to_fface_to_cface[d+1][fface] = cface
        d_to_fface_to_cface_dim[d+1][fface] = cface_dim
      end
    end
  end
end

# fine face_labeling

fine_tag_to_name = copy(coarse_labeling.tag_to_name)
fine_tag_to_entities = copy(coarse_labeling.tag_to_entities)

fine_d_to_dface_to_entity = Vector{Vector{Int32}}(undef,Dc+1)
for d in 0:Dc
  nF = num_faces(ftopo,d)
  fine_dface_to_entity = Vector{Int32}(undef,nF)

  for fface in 1:nF
    cface = d_to_fface_to_cface[d+1][fface]
    cface_dim = d_to_fface_to_cface_dim[d+1][fface]
    cface_entity = coarse_labeling.d_to_dface_to_entity[cface_dim+1][cface]
    fine_dface_to_entity[fface] = cface_entity
  end

  fine_d_to_dface_to_entity[d+1] = fine_dface_to_entity
end

fine_labeling = Geometry.FaceLabeling(fine_d_to_dface_to_entity,fine_tag_to_entities,fine_tag_to_name)

# Test implementation

vertex_coordinates = Geometry.get_node_coordinates(fmodel)
d_to_face_to_vertex = map(d->Geometry.get_faces(ftopo,d,0),0:Dc)

d_to_centroid = map(d_to_face_to_vertex) do face_to_vertex
  map(face_to_vertex) do vertices 
    sum(vertex_coordinates[vertices])/length(vertices)
  end
end

left_boundary_by_tags = findall(get_face_mask(fine_labeling,"left_boundary",1))
left_boundary_by_coords = findall(p->abs(p[1]) < 1.e-3,d_to_centroid[2])
@test sort(left_boundary_by_tags) == sort(left_boundary_by_coords)

bottom_corners_by_tags = findall(get_face_mask(fine_labeling,"bottom_corners",0))
bottom_corners_by_coords = findall(p->(abs(p[2]) < 1.e-3) && ((abs(p[1]) < 1.e-3) || (abs(p[1]-1.0) < 1.e-3)),d_to_centroid[1])
@test sort(bottom_corners_by_tags) == sort(bottom_corners_by_coords)

end