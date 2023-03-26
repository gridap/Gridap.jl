module FaceLabelingTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.Adaptivity

Dc = 2
model  = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))
ctopo = get_grid_topology(model)

coarse_labeling = get_face_labeling(model)
add_tag_from_tags!(coarse_labeling,"left_boundary",[7])
add_tag_from_tags!(coarse_labeling,"bottom_corners",[1,2])

fmodel = refine(model)
ftopo  = get_grid_topology(fmodel)
fine_labeling = get_face_labeling(fmodel)

# Compute face centroids 
vertex_coordinates = Geometry.get_node_coordinates(fmodel)
d_to_face_to_vertex = map(d->Geometry.get_faces(ftopo,d,0),0:Dc)

d_to_centroid = map(d_to_face_to_vertex) do face_to_vertex
  map(face_to_vertex) do vertices 
    sum(vertex_coordinates[vertices])/length(vertices)
  end
end

# Check tags
left_boundary_by_tags = findall(get_face_mask(fine_labeling,"left_boundary",1))
left_boundary_by_coords = findall(p->abs(p[1]) < 1.e-3,d_to_centroid[2])
@test sort(left_boundary_by_tags) == sort(left_boundary_by_coords)

bottom_corners_by_tags = findall(get_face_mask(fine_labeling,"bottom_corners",0))
bottom_corners_by_coords = findall(p->(abs(p[2]) < 1.e-3) && ((abs(p[1]) < 1.e-3) || (abs(p[1]-1.0) < 1.e-3)),d_to_centroid[1])
@test sort(bottom_corners_by_tags) == sort(bottom_corners_by_coords)



end