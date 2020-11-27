
struct DiracDelta{D} <: GridapType
  Γ::Triangulation{D}
  dΓ::LebesgueMeasure
end

function DiracDelta{D}(
  model::DiscreteModel,
  face_to_bgface::AbstractVector{<:Integer},
  degree::Integer) where D

  @assert 0 <= D && D < num_cell_dims(model) """\n
  Incorrect value of D=$D for building a DiracDelta{D} on a model with $(num_cell_dims(model)) cell dims.

  D should be in [0,$(num_cell_dims(model))).
  """
  topo = get_grid_topology(model)
  bgface_grid = Grid(ReferenceFE{D},model)
  face_trian = RestrictedTriangulation(bgface_grid,face_to_bgface)
  cell_trian = Triangulation(model)
  bgface_to_lcell = Fill(1,num_faces(model,D))
  glue = Geometry.FaceToCellGlue(topo,cell_trian,face_trian,face_to_bgface,bgface_to_lcell)
  Γ = BoundaryTriangulation(face_trian,cell_trian,glue)
  dΓ = LebesgueMeasure(Γ,degree)
  DiracDelta(Γ,dΓ)
end

function DiracDelta{D}(
  model::DiscreteModel,
  bgface_to_mask::AbstractVector{Bool},
  degree::Integer) where D

  face_to_bgface = findall(bgface_to_mask)
  DiracDelta{D}(model,face_to_bgface,degree)
end

function DiracDelta{D}(
  model::DiscreteModel,
  labeling::FaceLabeling,
  degree::Integer;tags) where D

  @assert 0 <= D && D < num_cell_dims(model) """\n
  Incorrect value of D=$D for building a DiracDelta{D} on a model with $(num_cell_dims(model)) cell dims.

  D should be in [0,$(num_cell_dims(model))).
  """
  face_to_mask = get_face_mask(labeling,tags,D)
  DiracDelta{D}(model,face_to_mask,degree)
end 

function DiracDelta{D}(model::DiscreteModel,degree::Integer;tags) where D
  labeling = get_face_labeling(model)
  DiracDelta{D}(model,labeling,degree,tags=tags)
end

function DiracDelta{0}(model::DiscreteModel,labeling::FaceLabeling;tags)
  degree = 0
  DiracDelta{0}(model,labeling,degree;tags=tags)
end

function DiracDelta{0}(model::DiscreteModel;tags)
  degree = 0
  DiracDelta{0}(model,degree;tags=tags)
end

function (d::DiracDelta)(f::CellField)
  evaluate(d,f)
end

function evaluate!(cache,d::DiracDelta,f::CellField)
  ∫(f)*d.dΓ
end

function get_triangulation(d::DiracDelta)
  d.Γ
end

