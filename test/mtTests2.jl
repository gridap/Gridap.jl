using Gridap
using Gridap.Geometry

D = 3
domain = repeat([0.0,1.0],D)
model = CartesianDiscreteModel(domain,(2,2,3))

face_tags = ["tag_23","tag_26","tag_24","tag_25"]
Γ  = BoundaryTriangulation(model,tags=face_tags)
_Λ  = SkeletonTriangulation(Γ)
glue = get_glue(_Λ,Val(D-1)) # glue from the edges to the attached faces

edge_tags = ["tag_17","tag_18","tag_20","tag_19"]
labeling = get_face_labeling(model)
edge_to_mask = get_face_mask(labeling,edge_tags,D-2)
edge_to_bgedge = findall(edge_to_mask) # ids of the edges we want

topo = get_grid_topology(model)
edge_to_face_map = Geometry.get_faces(topo,D-2,D-1) # edge to face connectivity

# Find local numbering of the edges we want
face_pairs = map(p -> [p...],zip(glue.plus.tface_to_mface,glue.minus.tface_to_mface))
local_edge_ids = map(edge_to_bgedge) do edge
  edge_faces = edge_to_face_map[edge]
  return findfirst(faces -> faces == edge_faces, face_pairs)
end
Λ = view(_Λ,local_edge_ids)

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(Γ,reffe)

dΛ = Measure(Λ,2)
v = interpolate(x -> sum(x),V)

l(v) = ∫(mean(v))*dΛ
b = assemble_vector(l,V)
