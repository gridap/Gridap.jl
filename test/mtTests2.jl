using Gridap
using Gridap.Geometry
using FillArrays

function EdgeTriangulation(model::DiscreteModel,tags)
  D = Geometry.num_cell_dims(model)
  labeling = get_face_labeling(model)
  face_to_mask = get_face_mask(labeling,tags,D-2)
  face_to_bgface = findall(face_to_mask)
  return EdgeTriangulation(model,face_to_bgface)
end

function EdgeTriangulation(
  model::DiscreteModel,
  face_to_bgface::AbstractVector{<:Integer})

  D = Geometry.num_cell_dims(model)
  topo = get_grid_topology(model)
  bgface_grid = Grid(Geometry.ReferenceFE{D-2},model)
  bgface_to_lcell = Fill(1,Geometry.num_faces(model,D-2))

  face_grid = view(bgface_grid,face_to_bgface)
  cell_grid = get_grid(model)
  glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
  trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

  BoundaryTriangulation(trian,glue)
end

model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,3))

face_tags = ["tag_23","tag_26","tag_24","tag_25"]
Γ  = BoundaryTriangulation(model,tags=face_tags)
nF = num_cells(Γ)

edge_tags = ["tag_17","tag_18","tag_20","tag_19"]
Λ  = EdgeTriangulation(model,edge_tags)
nE = num_cells(Λ)

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(Γ,reffe)

v = zero(V)
pts = get_cell_points(Λ)
v(pts)

