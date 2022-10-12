abstract type DiracDeltaSupportType <: GridapType end
struct IsGridEntity <: DiracDeltaSupportType end
struct NotGridEntity <: DiracDeltaSupportType end

const _is_grid_entity = IsGridEntity()
const _not_grid_entity = NotGridEntity()


struct DiracDelta{D,S<:DiracDeltaSupportType} <: GridapType
  Γ::Triangulation{D}
  dΓ::Measure
  supp::S
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
  face_grid = view(bgface_grid,face_to_bgface)
  cell_grid = get_grid(model)
  bgface_to_lcell = Fill(1,num_faces(model,D))
  glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
  trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)
  Γ = BoundaryTriangulation(trian,glue)
  dΓ = Measure(Γ,degree)
  DiracDelta(Γ,dΓ,_is_grid_entity)
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

function (d::DiracDelta)(f)
  evaluate(d,f)
end

function evaluate!(cache,d::DiracDelta,f)
  ∫(f)*d.dΓ
end

function get_triangulation(d::DiracDelta)
  d.Γ
end

# For handling DiracDelta at a generic Point in the domain #

function DiracDelta(x::Point{D}, model::DiscreteModel{D}) where D
  DiracDelta(x, model)
end

function DiracDelta(x::Point{D,T}, model::DiscreteModel{D}) where {D,T}
  # check if the point is inside an active cell, as it wouldn't be caught for
  # user-defined functions (i.e. which are not CellFields)
  trian = Triangulation(model)
  cache1 = _point_to_cell_cache(KDTreeSearch(),trian)
  cell = _point_to_cell!(cache1, x) # throws error if Point not in domain
  point_grid = UnstructuredGrid([x])
  point_model = UnstructuredDiscreteModel(point_grid)
  point_trian = Triangulation(point_model)
  dx = Measure(point_trian,1)
  DiracDelta(point_trian,dx,_not_grid_entity)
end

function DiracDelta(v::Vector{Point{D,T}},model::DiscreteModel{D}) where {D,T}
  # check if the points are inside an active cell
  trian = Triangulation(model)
  cache1 = _point_to_cell_cache(KDTreeSearch(),trian)
  cell = map(x -> _point_to_cell!(cache1, x), v)
  # throws error if any Point not in domain
  point_grid = UnstructuredGrid(v)
  point_model = UnstructuredDiscreteModel(point_grid)
  point_trian = Triangulation(point_model)
  dx = Measure(point_trian,1)
  DiracDelta(point_trian,dx,_not_grid_entity)
end

function evaluate!(cache,d::DiracDelta{0,NotGridEntity},f::CellField)
  eval = lazy_map(f,d.Γ.grid.node_coordinates)
  dc = DomainContribution()
  add_contribution!(dc, d.Γ, eval)
end
