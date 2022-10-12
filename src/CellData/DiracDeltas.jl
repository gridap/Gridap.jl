
struct DiracDelta{D} <: GridapType
  Γ::Triangulation{D}
  dΓ::Measure
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

#=
The idea would be to use the existing structure of DiracDelta.
In place of triangulation, we need something similar for Measure
it should be in PhysicalSpace to not repeat integration
reuse the cache from KDTree , if possible.

we also need to add the contribution to the respective cell, which apparently
was not done in the existing implementation of DiracDelta. One way to this is
to add the contribution to the respective cell where the point belongs to!
So that there is no problem with AD happening cell wise
=#
function DiracDelta(x::Point)
  DiracDelta{0}(x)
end

function DiracDelta{0}(x::Point{D,T}) where {D,T}
  desc = CartesianDescriptor(x, ntuple(i->zero(T),Val(D)), ntuple(i->1,Val(D)))
  point_model = CartesianDiscreteModel(desc)
  point_Ω = Triangulation(point_model)
  point_Γ = BoundaryTriangulation(point_model)
  quad = GenericQuadrature([x],ones(T,1),"DiracDelta point node")
  cell_quad = CellQuadrature(SA[quad],SA[quad.coordinates],SA[quad.weights],point_Ω, PhysicalDomain(),PhysicalDomain())
  dx = Measure(cell_quad)
  DiracDelta{0}(point_Γ,dx)
end

# goal is to build a DiscreteModel of just a single Point 
