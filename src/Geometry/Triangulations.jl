"""
    abstract type Triangulation{Dt,Dp}

A discredited physical domain associated with a `DiscreteModel{Dm,Dp}`.

`Dt` and `Dm` can be different.

The (mandatory) `Triangulation` interface can be tested with

- [`test_triangulation`](@ref)
"""
abstract type Triangulation{Dc,Dp} <: Grid{Dc,Dp} end

function get_discrete_model(t::Triangulation)
  @abstractmethod
end

function get_grid(t::Triangulation)
  @abstractmethod
end

# See possible types of glue below
function get_glue(t::Triangulation,::Val{d}) where d
  @abstractmethod
end

"""
    test_triangulation(trian::Triangulation)
"""
function test_triangulation(trian::Triangulation{Dc,Dp}) where {Dc,Dp}
  test_grid(trian)
  model = get_discrete_model(trian)
  @assert isa(model,DiscreteModel)
  grid = get_grid(trian)
  test_grid(grid)
end

# Grid interface
get_node_coordinates(trian::Triangulation) = get_node_coordinates(get_grid(trian))
get_cell_node_ids(trian::Triangulation) = get_cell_node_ids(get_grid(trian))
get_reffes(trian::Triangulation) = get_reffes(get_grid(trian))
get_cell_type(trian::Triangulation) = get_cell_type(get_grid(trian))
get_facet_normal(trian::Triangulation) = get_facet_normal(get_grid(trian))

# The following are not strictly needed, sine there is a default implementation for them.
# In any case, we delegate just in case the underlying grid defines more performant versions
get_cell_coordinates(trian::Triangulation) = get_cell_coordinates(get_grid(trian))
get_cell_shapefuns(trian::Triangulation) = get_cell_shapefuns(get_grid(trian))
get_cell_map(trian::Triangulation) = get_cell_map(get_grid(trian))
get_cell_reffe(trian::Triangulation) = get_cell_reffe(get_grid(trian))
is_first_order(trian::Triangulation) = is_first_order(get_grid(trian))

# This is the most used glue, but others are possible, see e.g. SkeletonGlue.
struct FaceToFaceGlue{A,B}
  tface_to_mface::A
  tface_to_mface_map::B
end

# This is the most basic Triangulation
# It represents a physical domain built using the faces of a DiscreteModel
struct BodyFittedTriangulation{Dt,Dp,A,B,C} <: Triangulation{Dt,Dp}
  model::A
  grid::B
  tface_to_mface::C
  function BodyFittedTriangulation(model::DiscreteModel,grid::Grid,tface_to_mface)
    Dp = num_point_dims(model)
    @assert Dp == num_point_dims(grid)
    Dt = num_cell_dims(grid)
    A = typeof(model)
    B = typeof(grid)
    C = typeof(tface_to_mface)
    new{Dt,Dp,A,B,C}(model,grid,tface_to_mface)
  end
end

get_discrete_model(trian::BodyFittedTriangulation) = trian.model
get_grid(trian::BodyFittedTriangulation) = trian.grid

function get_glue(trian::BodyFittedTriangulation{Dt},::Val{Dt}) where Dt
  tface_to_mface_map = Fill(GenericField(identity),num_cells(trian))
  FaceToFaceGlue(trian.tface_to_mface,tface_to_mface_map)
end

function get_glue(trian::BodyFittedTriangulation{Dt},::Val{Dm}) where {Dt,Dm}
  @notimplemented
end

get_triangulation(model) = Triangulation(model)

function Triangulation(
  ::Type{ReferenceFE{d}},model::DiscreteModel,filter::AbstractArray) where d
  mgrid = Grid(ReferenceFE{d},model)
  # Grid portion is OK here since this is usally used to
  # define a FE space
  tgrid = GridPortion(mgrid,filter)
  tface_to_mface = tgrid.cell_to_parent_cell
  BodyFittedTriangulation(model,tgrid,tface_to_mface)
end

function Triangulation(model::DiscreteModel,filter::AbstractArray)
  d = num_cell_dims(model)
  Triangulation(ReferenceFE{d},model,filter)
end

function Triangulation(
  ::Type{ReferenceFE{d}},
  model::DiscreteModel,
  labels::FaceLabeling;tags=nothing) where d

  if tags === nothing
    grid = Grid(ReferenceFE{d},model)
    tface_to_mface = IdentityVector(num_cells(grid))
    BodyFittedTriangulation(model,grid,tface_to_mface)
  else
    mface_to_mask = get_face_mask(labels,tags,d)
    Triangulation(ReferenceFE{d},model,mface_to_mask)
  end
end

function Triangulation(
  ::Type{ReferenceFE{d}},model::DiscreteModel;kwargs...) where d
  labels = get_face_labeling(model)
  Triangulation(ReferenceFE{d},model,labels;kwargs...)
end

function Triangulation(model::DiscreteModel;kwargs...)
  d = num_cell_dims(model)
  labels = get_face_labeling(model)
  Triangulation(ReferenceFE{d},model,labels;kwargs...)
end
