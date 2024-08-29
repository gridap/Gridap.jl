"""
    abstract type Triangulation{Dt,Dp}

A discredited physical domain associated with a `DiscreteModel{Dm,Dp}`.

`Dt` and `Dm` can be different.

The (mandatory) `Triangulation` interface can be tested with

- [`test_triangulation`](@ref)
"""
abstract type Triangulation{Dc,Dp} <: Grid{Dc,Dp} end

function get_background_model(t::Triangulation)
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
  model = get_background_model(trian)
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

"""
    is_change_possible(strian::Triangulation,ttrian::Triangulation)

  Returns `true` if `CellDatum` objects can be transferred from `strian` to `ttrian`.
"""
function is_change_possible(strian::Triangulation,ttrian::Triangulation)
  if strian === ttrian
    return true
  end
  @check get_background_model(strian) === get_background_model(ttrian) "Triangulations do not point to the same background discrete model!"
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

"""
    best_target(trian1::Triangulation,trian2::Triangulation)

  If possible, returns a `Triangulation` to which `CellDatum` objects can be transferred 
  from `trian1` and `trian2`. Can be `trian1`, `trian2` or a new `Triangulation`.
"""
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
  model = get_background_model(trian1)
  D = num_cell_dims(trian1)
  @assert num_cell_dims(trian2) == D
  Triangulation(ReferenceFE{D},model)
end

function get_active_model(t::Triangulation)
  compute_active_model(t)
end

function compute_active_model(t::Triangulation)
  D = num_cell_dims(t)
  glue = get_glue(t,Val(D))
  @assert glue.mface_to_tface !== nothing
  bgmodel = get_background_model(t)
  model = DiscreteModel(Polytope{D},bgmodel)
  _restrict_model(model,get_grid(t),glue.tface_to_mface)
end

function _restrict_model(model,grid::Grid,tface_to_mface)
  _restrict_model(model,tface_to_mface)
end

function _restrict_model(model,grid::GridPortion,tface_to_mface)
  @check grid.cell_to_parent_cell == tface_to_mface
  DiscreteModelPortion(model,grid)
end

function _restrict_model(model,tface_to_mface)
  DiscreteModelPortion(model,tface_to_mface)
end

function _restrict_model(model,tface_to_mface::IdentityVector)
  model
end

abstract type TrianFaceModelFaceMapInjectivity end;
struct Injective    <: TrianFaceModelFaceMapInjectivity end;
struct NonInjective <: TrianFaceModelFaceMapInjectivity end;

# This is the most basic Triangulation
# It represents a physical domain built using the faces of a DiscreteModel
struct BodyFittedTriangulation{Dt,Dp,A,B,C,D<:TrianFaceModelFaceMapInjectivity} <: Triangulation{Dt,Dp}
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

    # While we do not have a more definitive solution, we need to distinguish 
    # between injective and non-injective tface_to_mface maps.
    # The inverse map, mface_to_tface, relies on PosNegPartition, which fails 
    # whenever the same mface is the image of more than one tface.
    # In turn, I have required non-injective mappings for the computation of facet 
    # integrals on non-conforming cell interfaces.
    if !(allunique(tface_to_mface))
      tface_to_mface_injectivity = NonInjective()
      D = typeof(tface_to_mface_injectivity)
      new{Dt,Dp,A,B,C,D}(model,grid,tface_to_mface)
    else 
      tface_to_mface_injectivity = Injective()
      D = typeof(tface_to_mface_injectivity)
      new{Dt,Dp,A,B,C,D}(model,grid,tface_to_mface)
    end
  end
end

get_background_model(trian::BodyFittedTriangulation) = trian.model
get_grid(trian::BodyFittedTriangulation) = trian.grid

function get_glue(trian::BodyFittedTriangulation{Dt,Dp,A,B,C,Injective},::Val{Dt}) where {Dt,Dp,A,B,C}
  tface_to_mface_map = Fill(GenericField(identity),num_cells(trian))
  if isa(trian.tface_to_mface,IdentityVector) && num_faces(trian.model,Dt) == num_cells(trian)
    mface_to_tface = trian.tface_to_mface
  else
    nmfaces = num_faces(trian.model,Dt)
    mface_to_tface = PosNegPartition(trian.tface_to_mface,Int32(nmfaces))
  end
  FaceToFaceGlue(trian.tface_to_mface,tface_to_mface_map,mface_to_tface)
end

function get_glue(trian::BodyFittedTriangulation{Dt,Dp,A,B,C,NonInjective},::Val{Dt}) where {Dt,Dp,A,B,C}
  tface_to_mface_map = Fill(GenericField(identity),num_cells(trian))
  mface_to_tface = nothing
  # Whenever tface_to_mface is non-injective, we currently avoid the computation of 
  # mface_to_tface, which relies on PosNegPartition. This is a limitation that we should 
  # face in the future on those scenarios on which we need mface_to_tface.
  FaceToFaceGlue(trian.tface_to_mface,tface_to_mface_map,mface_to_tface)
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
  # Grid portion is OK here since this is usually used to
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

function Triangulation(trian::Triangulation,args...;kwargs...)
  amodel = get_active_model(trian)
  dtrian = Triangulation(amodel,args...;kwargs...)
  CompositeTriangulation(trian,dtrian)
end

function Triangulation(trian::Triangulation)
  trian
end

function Triangulation(trian::Triangulation,x::AbstractArray{<:Integer})
  view(trian,x)
end

function Triangulation(trian::Triangulation,x::AbstractArray{<:Bool})
  y = findall(collect1d(x))
  view(trian,y)
end

function Interior(args...;kwargs...)
  Triangulation(args...;kwargs...)
end

# This is the low-level functionality to move from one Triangulation to another

function restrict(a::AbstractArray,b::AbstractArray)
  lazy_map(Reindex(a),b)
end

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

function extend(a::LazyArray{<:Fill{typeof(transpose)}},b::PosNegPartition)
  c = a.args[1]
  d = extend(c,b)
  lazy_map(transpose,d)
end

function extend(a::LazyArray{<:Fill{typeof(linear_combination)}},b::PosNegPartition)
  d1 = extend(a.args[1],b)
  d2 = extend(a.args[2],b)
  lazy_map(linear_combination,d1,d2)
end

#function extend(a::LazyArray{<:Fill},b::PosNegPartition)
#  k = a.maps.value
#  args = map(i->extend(i,b),a.args)
#  lazy_map(k,args...)
#end

"""
"""
function pos_neg_data(ipos_to_val::AbstractArray,i_to_iposneg::PosNegPartition)
  @abstractmethod
end

function pos_neg_data(
  ipos_to_val::AbstractArray{<:Number},i_to_iposneg::PosNegPartition)
  nineg = length(i_to_iposneg.ineg_to_i)
  ineg_to_val = Fill(zero(eltype(ipos_to_val)),nineg)
  ipos_to_val, ineg_to_val
end

function pos_neg_data(
  ipos_to_val::AbstractArray{<:AbstractArray{<:Number}},i_to_iposneg::PosNegPartition)
  nineg = length(i_to_iposneg.ineg_to_i)
  val = testitem(ipos_to_val)
  zs = 0 .* size(val)
  void = similar(val,eltype(val),zs)
  ineg_to_val = Fill(void,nineg)
  ipos_to_val, ineg_to_val
end

function pos_neg_data(
  ipos_to_val::AbstractArray{<:Field},i_to_iposneg::PosNegPartition)
  nineg = length(i_to_iposneg.ineg_to_i)
  ipos_to_v = lazy_map(VoidFieldMap(false),ipos_to_val)
  ineg_to_v = Fill(VoidField(testitem(ipos_to_val),true),nineg)
  ipos_to_v, ineg_to_v
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

function pos_neg_data(
  ipos_to_val::AbstractArray{<:ArrayBlock},i_to_iposneg::PosNegPartition)
  nineg = length(i_to_iposneg.ineg_to_i)
  val = testitem(ipos_to_val)
  void = _similar_empty(val)
  ineg_to_val = Fill(void,nineg)
  ipos_to_val, ineg_to_val
end

function pos_neg_data(
  ipos_to_val::AbstractArray{<:Tuple{<:Any,<:Any}},i_to_iposneg::PosNegPartition)
  nineg = length(i_to_iposneg.ineg_to_i)
  val = testitem(ipos_to_val)
  void = _similar_empty(val)
  ineg_to_val = Fill(void,nineg)
  ipos_to_val, ineg_to_val
end

function _similar_empty(val::AbstractArray)
  zs = 0 .* size(val)
  void = similar(val,eltype(val),zs)
end

function _similar_empty(val::ArrayBlock)
  a = deepcopy(val)
  for i in eachindex(a)
    if a.touched[i]
      a.array[i] = _similar_empty(a.array[i])
    end
  end
  a
end

function _similar_empty(val::Tuple)
  a, b = val
  a1 = _similar_empty(a)
  b1 = _similar_empty(b)
  (a1,b1)
end


# "Compose" triangulations

struct CompositeTriangulation{Dc,Dp,A,B} <: Triangulation{Dc,Dp}
  rtrian::A
  dtrian::B
  function CompositeTriangulation(
    rtrian::Triangulation, dtrian::Triangulation)
    @assert num_point_dims(rtrian) == num_point_dims(dtrian)
    Dp = num_point_dims(rtrian)
    Dc = num_cell_dims(dtrian)
    #@assert get_active_model(rtrian) === get_background_model(dtrian)
    A = typeof(rtrian)
    B = typeof(dtrian)
    new{Dc,Dp,A,B}(rtrian,dtrian)
  end
end

get_background_model(t::CompositeTriangulation) = get_background_model(t.rtrian)
get_active_model(t::CompositeTriangulation) = get_active_model(t.dtrian)
get_grid(t::CompositeTriangulation) = get_grid(t.dtrian)
get_facet_normal(t::CompositeTriangulation) = get_facet_normal(t.dtrian)
function get_glue(t::CompositeTriangulation,::Val{D}) where D
  Dr = num_cell_dims(t.rtrian)
  rglue = get_glue(t.rtrian,Val(D))
  dglue = get_glue(t.dtrian,Val(Dr))
  _compose_glues(rglue,dglue)
end

function _compose_glues(rglue,dglue)
  nothing
end

function _compose_glues(rglue::FaceToFaceGlue,dglue::FaceToFaceGlue)
  rface_to_mface = rglue.tface_to_mface
  dface_to_rface = dglue.tface_to_mface
  dface_to_mface = collect(lazy_map(Reindex(rface_to_mface),dface_to_rface))
  rface_to_mface_map = rglue.tface_to_mface_map
  dface_to_rface_map = dglue.tface_to_mface_map
  dface_to_mface_map1 = lazy_map(Reindex(rface_to_mface_map),dface_to_rface)
  dface_to_mface_map = lazy_map(∘,dface_to_mface_map1,dface_to_rface_map)
  mface_to_dface = nothing
  FaceToFaceGlue(dface_to_mface,dface_to_mface_map,mface_to_dface)
end

struct GenericTriangulation{Dc,Dp,A,B,C} <: Triangulation{Dc,Dp}
  grid::A
  model::B
  glue::C
  function GenericTriangulation(
    grid::Grid,
    model=nothing,
    glue=ntuple(i->nothing,num_cell_dims(grid)+1))
    Dc = num_cell_dims(grid)
    Dp = num_point_dims(grid)
    A = typeof(grid)
    B = typeof(model)
    C = typeof(glue)
    new{Dc,Dp,A,B,C}(grid,model,glue)
  end
end

get_grid(a::GenericTriangulation) = a.grid
get_glue(a::GenericTriangulation,::Val{D}) where D = a.glue[D+1]
function get_background_model(a::GenericTriangulation)
  @notimplementedif a.model === nothing "This triangulation object cannot be used to define a FE Space"
  a.model
end

struct TriangulationView{Dc,Dp,A,B} <: Triangulation{Dc,Dp}
  parent::A
  cell_to_parent_cell::B
  function TriangulationView(parent::Triangulation,cell_to_parent_cell::AbstractArray)
    Dc = num_cell_dims(parent)
    Dp = num_point_dims(parent)
    A = typeof(parent)
    B = typeof(cell_to_parent_cell)
    new{Dc,Dp,A,B}(parent,cell_to_parent_cell)
  end
end

Base.view(a::Triangulation,b::AbstractArray) = TriangulationView(a,b)

function TriangulationView(parent::Triangulation,parent_cell_to_mask::AbstractArray{Bool})
  cell_to_parent_cell = findall(collect1d(parent_cell_to_mask))
  TriangulationView(parent,cell_to_parent_cell)
end

function TriangulationView(parent::Triangulation,parent_cell_to_mask::AbstractVector{Bool})
  cell_to_parent_cell = findall(parent_cell_to_mask)
  TriangulationView(parent,cell_to_parent_cell)
end

function get_background_model(trian::TriangulationView)
  get_background_model(trian.parent)
end

function get_grid(trian::TriangulationView)
  view(get_grid(trian.parent),trian.cell_to_parent_cell)
end

function get_glue(t::TriangulationView,::Val{d}) where d
  parent = get_glue(t.parent,Val(d))
  parent === nothing && return nothing
  view(parent,t.cell_to_parent_cell)
end

function Base.view(glue::FaceToFaceGlue,ids::AbstractArray)
  tface_to_mface = lazy_map(Reindex(glue.tface_to_mface),ids)
  tface_to_mface_map = lazy_map(Reindex(glue.tface_to_mface_map),ids)
  if glue.mface_to_tface === nothing
    mface_to_tface = nothing
  else
    nmfaces = length(glue.mface_to_tface)
    mface_to_tface = PosNegPartition(tface_to_mface,Int32(nmfaces))
  end
  FaceToFaceGlue(tface_to_mface,tface_to_mface_map,mface_to_tface)
end

function get_facet_normal(trian::TriangulationView)
  n = get_facet_normal(trian.parent)
  restrict(n,trian.cell_to_parent_cell)
end

function get_cell_map(trian::TriangulationView)
  lazy_map(Reindex(get_cell_map(trian.parent)),trian.cell_to_parent_cell)
end
