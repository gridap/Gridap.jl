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
  nothing
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
get_cell_ref_coordinates(trian::Triangulation) = get_cell_ref_coordinates(get_grid(trian))
get_cell_shapefuns(trian::Triangulation) = get_cell_shapefuns(get_grid(trian))
get_cell_map(trian::Triangulation) = get_cell_map(get_grid(trian))
get_cell_reffe(trian::Triangulation) = get_cell_reffe(get_grid(trian))
is_first_order(trian::Triangulation) = is_first_order(get_grid(trian))

# This is the most used glue, but others are possible, see e.g. SkeletonGlue.
struct FaceToFaceGlue{A,B,C}
  tface_to_mface::A
  tface_to_mface_map::B
  mface_to_tface::C
end

function is_change_possible(strian::Triangulation,ttrian::Triangulation)
  @check strian !== ttrian "This function is not meant to be called on the same object"
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  is_change_possible(sglue,tglue)
end

function is_change_possible(strian,ttrian)
  false
end

function is_change_possible(sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  sglue.mface_to_tface != nothing
end

function best_target(trian1::Triangulation,trian2::Triangulation)
  @check is_change_possible(trian1,trian2)
  @check is_change_possible(trian2,trian1)
  D1 = num_cell_dims(trian1)
  D2 = num_cell_dims(trian2)
  glue1 = get_glue(trian1,Val(D2))
  glue2 = get_glue(trian2,Val(D1))
  best_target(trian1,trian2,glue1,glue2)
end

function best_target(
  trian1::Triangulation,trian2::Triangulation,glue1::FaceToFaceGlue,glue2::FaceToFaceGlue)
  # Return the background
  model = get_discrete_model(trian1)
  D = num_cell_dims(trian1)
  @assert num_cell_dims(trian2) == D
  Triangulation(ReferenceFE{D},model)
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
  tface_to_mface = trian.tface_to_mface
  tface_to_mface_map = Fill(GenericField(identity),num_cells(trian))
  if isa(tface_to_mface,IdentityVector) && num_faces(trian.model,Dt) == num_cells(trian)
    mface_to_tface = tface_to_mface
  else
    ntfaces = num_faces(trian.model,Dt)
    mface_to_tface = PosNegPartition(tface_to_mface,Int32(ntfaces))
  end
  FaceToFaceGlue(tface_to_mface,tface_to_mface_map,mface_to_tface)
end

#function get_glue(trian::BodyFittedTriangulation{Dt},::Val{Dm}) where {Dt,Dm}
#  @notimplemented
#end

function Base.view(t::BodyFittedTriangulation,ids::AbstractArray)
  model = t.model
  grid = view(t.grid,ids)
  tface_to_mface = lazy_map(Reindex(t.tface_to_mface),ids)
  BodyFittedTriangulation(model,grid,tface_to_mface)
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

# This is the low-level functionallity to move from one Triangulation to another

function extend(tface_to_val,mface_to_tface)
  @notimplemented
end

function extend(tface_to_val,mface_to_tface::IdentityVector)
  tface_to_val
end

function extend(tface_to_val,mface_to_tface::PosNegPartition)
  ipos_to_val, ineg_to_val = pos_neg_data(tface_to_val,mface_to_tface)
  i_to_iposneg = mface_to_tface
  lazy_map(PosNegReindex(ipos_to_val,ineg_to_val),i_to_iposneg)
end

function pos_neg_data(
  ipos_to_val::AbstractArray{<:Number},i_to_iposneg::PosNegPartition)
  nineg = length(i_to_iposneg.ineg_to_i)
  ineg_to_val = Fill(zero(eltype(ipos_to_val)),nineg)
  ipos_to_val, ineg_to_val
end

function pos_neg_data(
  ipos_to_val::AbstractArray{<:AbstractArray{<:Field}},i_to_iposneg::PosNegPartition)
  _pos_neg_data_basis(ipos_to_val,i_to_iposneg)
end

function pos_neg_data(
  ipos_to_val::AbstractArray{<:AbstractArray{<:Dof}},i_to_iposneg::PosNegPartition)
  _pos_neg_data_basis(ipos_to_val,i_to_iposneg)
end

function _pos_neg_data_basis(ipos_to_val,i_to_iposneg)
  nineg = length(i_to_iposneg.ineg_to_i)
  ipos_to_v = lazy_map(VoidBasisMap(false),ipos_to_val)
  ineg_to_v = Fill(VoidBasis(testitem(ipos_to_val),true),nineg)
  ipos_to_v, ineg_to_v
end

