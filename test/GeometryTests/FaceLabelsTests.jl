module FaceLabelsTests

using Test
using Gridap

vertex_to_geolabel = [1,1,2,2,2,1,1,3,3]
edge_to_geolabel = [4,4,5,5,5,5,6,6,4]
physlabel_1 = [1,3,4]
physlabel_2 = [5,3,6,2]
tag_to_name = ["label1","label2"]

labels = FaceLabels(
  [vertex_to_geolabel, edge_to_geolabel],
  [physlabel_1, physlabel_2],
  tag_to_name)

@test isa(labels,FaceLabels)
@test labels_on_dim(labels,0) == vertex_to_geolabel
@test labels_on_dim(labels,1) == edge_to_geolabel
@test labels_on_tag(labels,1) == physlabel_1
@test labels_on_tag(labels,2) == physlabel_2
@test tag_from_name(labels,"label1")==1

r = [1, 1, 2, 2, 2, 2, 2, 2, 1]
@test collect(first_tag_on_face(labels,1)) == r

end # module
