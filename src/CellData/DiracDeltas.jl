abstract type DiracDeltaSupportType <: GridapType end
abstract type IsGridEntity <: DiracDeltaSupportType end
abstract type NotGridEntity <: DiracDeltaSupportType end

struct GenericDiracDelta{D,Dt,S<:DiracDeltaSupportType} <: GridapType
  Γ::Triangulation{Dt}
  dΓ::Measure
end

"""
    const DiracDelta{D} = GenericDiracDelta{D,D,IsGridEntity}
"""
const DiracDelta{D} = GenericDiracDelta{D,D,IsGridEntity}

"""
    DiracDelta{D}(model::DiscreteModel, degree; tags)
    DiracDelta{0}(model::DiscreteModel; tags)
    DiracDelta{D}(model::DiscreteModel, labeling::FaceLabeling, degree; tags)
    DiracDelta{0}(model::DiscreteModel, labeling::FaceLabeling; tags)
    DiracDelta{D}(model::DiscreteModel, face_to_bgface::AbstractVector{<:Integer}, degree)
    DiracDelta{D}(model::DiscreteModel, bgface_to_mask::AbstractVector{Bool}, degree)
    DiracDelta(   model::DiscreteModel{D}, p::Point{D,T})
    DiracDelta(   model::DiscreteModel{D}, pvec::Vector{Point{D,T}})

where `degree` isa `Integer`.
"""
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
  DiracDelta{D}(Γ,dΓ)
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

function (d::GenericDiracDelta)(f)
  evaluate(d,f)
end

function evaluate!(cache,d::GenericDiracDelta,f)
  ∫(f)*d.dΓ
end

function get_triangulation(d::GenericDiracDelta)
  d.Γ
end

# For handling DiracDelta at a generic Point in the domain #

function _cell_to_pindices(pvec::Vector{<:Point},trian::Triangulation)
  cache = _point_to_cell_cache(KDTreeSearch(),trian)
  cell_to_pindex = Dict{Int, Vector{Int32}}()
  for i in 1:length(pvec)
    cell = _point_to_cell!(cache, pvec[i])
    push!(get!(() -> valtype(cell_to_pindex)[], cell_to_pindex, cell), i)
  end
  cell_to_pindex
end

function DiracDelta(model::DiscreteModel{D}, p::Point{D,T}) where {D,T}
  trian = Triangulation(model)
  cache = _point_to_cell_cache(KDTreeSearch(),trian)
  cell = _point_to_cell!(cache, p)
  trianv = view(trian,[cell])
  point = [p]
  weight = [one(T)]
  pquad = GenericQuadrature(point,weight)
  pmeas = Measure(CellQuadrature([pquad],[pquad.coordinates],[pquad.weights],trianv,PhysicalDomain(),PhysicalDomain()))
  GenericDiracDelta{0,D,NotGridEntity}(trianv,pmeas)
end

function DiracDelta(model::DiscreteModel{D}, pvec::Vector{Point{D,T}}) where {D,T}
  trian = Triangulation(model)
  cell_to_pindices = _cell_to_pindices(pvec,trian)
  cell_ids = collect(keys(cell_to_pindices))
  cell_points = collect(values(cell_to_pindices))
  points = map(i->pvec[cell_points[i]], 1:length(cell_ids))
  weights_x_cell = Fill.(one(T),length.(cell_points))
  pquad = map(i -> GenericQuadrature(points[i],weights_x_cell[i]), 1:length(cell_ids))
  trianv = view(trian,cell_ids)
  pmeas = Measure(CellQuadrature(pquad,points,weights_x_cell,trianv,PhysicalDomain(),PhysicalDomain()))
  GenericDiracDelta{0,D,NotGridEntity}(trianv,pmeas)
end
