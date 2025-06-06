"""
A single point or an array of points on the cells of a Triangulation
CellField objects can be evaluated efficiently at CellPoint instances.
"""
struct CellPoint{DS,A,B,C} <: CellDatum
  cell_ref_point::A
  cell_phys_point::B
  trian::C
  domain_style::DS
end

function CellPoint(
  cell_ref_point::AbstractArray,
  trian::Triangulation,
  domain_style::ReferenceDomain)

  cell_map = get_cell_map(trian)
  cell_phys_point = lazy_map(evaluate,cell_map,cell_ref_point)
  CellPoint(cell_ref_point,cell_phys_point,trian,domain_style)
end

function CellPoint(
  cell_phys_point::AbstractArray,
  trian::Triangulation,
  domain_style::PhysicalDomain)
  cell_map = get_cell_map(trian)
  cell_invmap = lazy_map(inverse_map,cell_map)
  cell_ref_point = lazy_map(evaluate,cell_invmap,cell_phys_point)
  CellPoint(cell_ref_point,cell_phys_point,trian,domain_style)
end

function get_data(f::CellPoint)
  if DomainStyle(f) == ReferenceDomain()
    f.cell_ref_point
  else
    f.cell_phys_point
  end
end

get_triangulation(f::CellPoint) = f.trian
DomainStyle(::Type{CellPoint{DS,A,B,C}}) where {DS,A,B,C} = DS()

function change_domain(a::CellPoint,::ReferenceDomain,::PhysicalDomain)
  CellPoint(a.cell_ref_point,a.cell_phys_point,a.trian,PhysicalDomain())
end

function change_domain(a::CellPoint,::PhysicalDomain,::ReferenceDomain)
  CellPoint(a.cell_ref_point,a.cell_phys_point,a.trian,ReferenceDomain())
end

# Possibly with a different name
"""
"""
function get_cell_points(trian::Triangulation)
  cell_ref_coords = get_cell_ref_coordinates(trian)
  cell_phys_coords = get_cell_coordinates(trian)
  CellPoint(cell_ref_coords,cell_phys_coords,trian,ReferenceDomain())
end

function Base.:(==)(a::CellPoint,b::CellPoint)
  a.trian == b.trian &&
  a.cell_ref_point == b.cell_ref_point &&
  a.cell_phys_point == b.cell_phys_point &&
  a.domain_style == b.domain_style
end

"""
"""
abstract type CellField <: CellDatum end

function similar_cell_field(f::CellField,cell_data,trian,ds)
  GenericCellField(cell_data,trian,ds)
end

function Base.show(io::IO,::MIME"text/plain",f::CellField)
  show(io,f)
  print(io,":")
  print(io,"\n num_cells: $(num_cells(f))")
  print(io,"\n DomainStyle: $(DomainStyle(f))")
  print(io,"\n Triangulation: $(get_triangulation(f))")
  print(io,"\n Triangulation id: $(objectid(get_triangulation(f)))")
end

function CellField(f::Function,trian::Triangulation,domain_style::DomainStyle)
  s = size(get_cell_map(trian))
  cell_field = Fill(GenericField(f),s)
  GenericCellField(cell_field,trian,PhysicalDomain())
end

function CellField(f::Number,trian::Triangulation,domain_style::DomainStyle)
  s = size(get_cell_map(trian))
  cell_field = Fill(constant_field(f),s)
  GenericCellField(cell_field,trian,domain_style)
end

function CellField(f::AbstractArray{<:Number},trian::Triangulation,domain_style::DomainStyle)
  @check length(f)==num_cells(trian)  """\n
  You are trying to build a CellField from an array of length $(length(f))
  on a Triangulation with $(num_cells(trian)) cells. The length of the given array
  and the number of cells should match.
  """
  cell_field = lazy_map(constant_field,f)
  GenericCellField(cell_field,trian,domain_style)
end

function CellField(f::CellField,trian::Triangulation,domain_style::DomainStyle)
  change_domain(f,trian,domain_style)
end

function CellField(f,trian::Triangulation)
  CellField(f,trian,ReferenceDomain())
end

"""
    get_normal_vector(trian::Triangulation)
"""
function get_normal_vector(trian::Triangulation)
  cell_normal = get_facet_normal(trian)
  get_normal_vector(trian, cell_normal)
end

function get_tangent_vector(trian::Triangulation)
  cell_tangent = get_edge_tangent(trian)
  get_tangent_vector(trian, cell_tangent)
end

function get_normal_vector(trian::Triangulation,cell_vectors::AbstractArray)
  GenericCellField(cell_vectors,trian,ReferenceDomain())
end

function get_tangent_vector(trian::Triangulation,cell_vectors::AbstractArray)
  GenericCellField(cell_vectors,trian,ReferenceDomain())
end

evaluate!(cache,f::Function,x::CellPoint) = CellField(f,get_triangulation(x))(x)

function change_domain(a::CellField,::ReferenceDomain,::PhysicalDomain)
  trian = get_triangulation(a)
  cell_map = get_cell_map(trian)
  cell_invmap = lazy_map(inverse_map,cell_map)
  cell_field_ref = get_data(a)
  cell_field_phys = lazy_map(Broadcasting(∘),cell_field_ref,cell_invmap)
  similar_cell_field(a,cell_field_phys,trian,PhysicalDomain())
end

function change_domain(a::CellField,::PhysicalDomain,::ReferenceDomain)
  trian = get_triangulation(a)
  cell_map = get_cell_map(trian)
  cell_field_phys = get_data(a)
  cell_field_ref = lazy_map(Broadcasting(∘),cell_field_phys,cell_map)
  similar_cell_field(a,cell_field_ref,trian,ReferenceDomain())
end

function change_domain(a::CellField,target_trian::Triangulation,target_domain::DomainStyle)
  change_domain(a,get_triangulation(a),DomainStyle(a),target_trian,target_domain)
end

function change_domain(a::CellField,source_domain::DomainStyle,target_trian::Triangulation,target_domain::DomainStyle)
  change_domain(a,get_triangulation(a),source_domain,target_trian,target_domain)
end

function change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,ttrian::Triangulation,::ReferenceDomain)
  msg = """\n
  We cannot move the given CellField to the reference domain of the requested triangulation.
  Make sure that the given triangulation is either the same as the triangulation on which the
  CellField is defined, or that the latter triangulation is the background of the former.
  """
  if strian === ttrian
    return a
  end
  @check is_change_possible(strian,ttrian) msg
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  change_domain_ref_ref(a,ttrian,sglue,tglue)
end

function change_domain(a::CellField,strian::Triangulation,::PhysicalDomain,ttrian::Triangulation,::PhysicalDomain)
  msg = """\n
  We cannot move the given CellField to the physical domain of the requested triangulation.
  Make sure that the given triangulation is either the same as the triangulation on which the
  CellField is defined, or that the latter triangulation is the background of the former.
  """
  if strian === ttrian
    return a
  end
  @check is_change_possible(strian,ttrian) msg
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  change_domain_phys_phys(a,ttrian,sglue,tglue)
end

function change_domain(a::CellField,strian::Triangulation,::PhysicalDomain,trian::Triangulation,::ReferenceDomain)
  a_trian = change_domain(a,trian,PhysicalDomain())
  change_domain(a_trian,ReferenceDomain())
end

function change_domain(a::CellField,strian::Triangulation,::ReferenceDomain,trian::Triangulation,::PhysicalDomain)
  a_phys = change_domain(a,PhysicalDomain())
  change_domain(a_phys,trian,PhysicalDomain())
end

function change_domain_ref_ref(
  a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  sface_to_field = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  tface_to_mface_map = tglue.tface_to_mface_map
  mface_to_field = extend(sface_to_field,mface_to_sface)
  tface_to_field_s = lazy_map(Reindex(mface_to_field),tface_to_mface)
  tface_to_field_t = lazy_map(Broadcasting(∘),tface_to_field_s,tface_to_mface_map)
  similar_cell_field(a,tface_to_field_t,ttrian,ReferenceDomain())
end

function change_domain_phys_phys(
  a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  sface_to_field = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  mface_to_field = extend(sface_to_field,mface_to_sface)
  tface_to_field = lazy_map(Reindex(mface_to_field),tface_to_mface)
  similar_cell_field(a,tface_to_field,ttrian,PhysicalDomain())
end

"""
"""
struct GenericCellField{DS} <: CellField
  cell_field::AbstractArray
  trian::Triangulation
  domain_style::DS
  function GenericCellField(
    cell_field::AbstractArray,
    trian::Triangulation,
    domain_style::DomainStyle)

    DS = typeof(domain_style)
    new{DS}(Fields.MemoArray(cell_field),trian,domain_style)
  end
end

get_data(f::GenericCellField) = f.cell_field
get_triangulation(f::GenericCellField) = f.trian
DomainStyle(::Type{GenericCellField{DS}}) where DS = DS()
function similar_cell_field(f::GenericCellField,cell_data,trian,ds)
  GenericCellField(cell_data,trian,ds)
end

"""
   dist = distance(polytope::ExtrusionPolytope,
                   inv_cmap::Field,
                   x::Point)

Calculate distance from point `x` to the polytope. The polytope is
given by its type and by the inverse cell map, i.e. by the map from
the physical to the reference space.

Positive distances are outside the polytope, negative distances are
inside the polytope.

The distance is measured in an unspecified norm, currently the L∞
norm.
"""
function distance(polytope::ExtrusionPolytope, inv_cmap::Field, x::Point)
  extrusion = polytope.extrusion
  isempty(extrusion) && return zero(eltype(x))
  p = inv_cmap(x)
  if all(e == HEX_AXIS for e in extrusion)
    # Boundaries are at `a=0` and `a=1` in each direction
    return maximum(max(0 - a, a - 1) for a in p)
  elseif all(e == TET_AXIS for e in extrusion)
    # Calculate barycentric coordinates
    λ = Point(p..., 1 - sum(p))
    return maximum(-λ)
  else
    @notimplemented "Only hypercubes and simplices are implemented so far"
  end
end

function _point_to_cell!(cache, x::Point)
  searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache

  # Loop over the first m.num_nearest_vertex
  for (id,dist) in zip(knn(kdtree, SVector(Tuple(x)), searchmethod.num_nearest_vertices, true)...)

    # Find all neighbouring cells
    cells = getindex!(table_cache,vertex_to_cells,id)
    @assert !isempty(cells)

    # Calculate the distance from the point to all the cells. Without
    # round-off, and with non-overlapping cells, the distance would be
    # negative for exactly one cell and positive for all other ones. Due
    # to round-off, the distance can be slightly negative or slightly
    # positive for points near cell boundaries, in particular near
    # vertices. In this case, choose the cell with the smallest
    # distance, and check that the distance (if positive) is at most at
    # round-off level.
    T = eltype(dist)
    function cell_distance(cell::Integer)
      ctype = cell_to_ctype[cell]
      polytope = ctype_to_polytope[ctype]
      cmap = cell_map[cell]
      inv_cmap = inverse_map(cmap)
      return distance(polytope, inv_cmap, x)
    end
    # findmin, without allocating an array
    cell = zero(eltype(cells))
    dist = T(Inf)
    for jcell in cells
      jdist = cell_distance(jcell)
      if jdist < dist
        cell = jcell
        dist = jdist
      end
    end

    dist ≤ 1000eps(T) && return cell

  end

  # Output error message if cell not found
  @check false "Point $x is not inside any active cell"
end

function evaluate!(cache,f::CellField,x::Point)
  cache1,cache2 = cache
  cell_f_cache, f_cache, cell_f, f₀ = cache2
  @check f === f₀ "Wrong cache"

  cell = _point_to_cell!(cache1, x)
  cf = getindex!(cell_f_cache, cell_f, cell)
  fx = evaluate!(f_cache, cf, x)
  return fx
end

return_cache(f::CellField,xs::AbstractVector{<:Point}) = return_cache(f,testitem(xs))

# # Simple version:
# function evaluate!(cache,f::CellField,xs::AbstractVector{<:Point})
#   return map(x->evaluate!(cache,f,x), xs)
# end

"""
    make_inverse_table(i2j::AbstractVector{<:Integer}, nj::Int)
"""
function make_inverse_table(i2j::AbstractVector{<:Integer},nj::Int)
  ni = length(i2j)
  @assert nj≥0

  p = sortperm(i2j)
  # i2j[p] is a sorted array of js

  data = p
  ptrs = Array{Int}(undef, nj+1)
  i2lis = Array{Int}(undef, ni)
  jprev = zero(eltype(i2j))
  liprev = 0
  for (n,i) in enumerate(p)
    j = i2j[i]
    @assert jprev≤j≤nj
    li = (j==jprev ? liprev : 0) + 1
    ptrs[jprev+1:j] .= n
    i2lis[i] = li
    jprev = j
    liprev = li
  end
  ptrs[jprev+1:nj+1] .= ni+1
  j2is = Table(data,ptrs)

  return j2is,i2lis
end

# Efficient version:
function evaluate!(cache,f::CellField,point_to_x::AbstractVector{<:Point})
  cache1,cache2 = cache
  searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map = cache1
  cell_f_cache, f_cache, cell_f, f₀ = cache2
  @check f === f₀ "Wrong cache"

  ncells = length(cell_map)
  x_to_cell(x) = _point_to_cell!(cache1,x)
  point_to_cell = map(x_to_cell,point_to_x)
  cell_to_points,point_to_lpoint = make_inverse_table(point_to_cell,ncells)
  cell_to_xs = lazy_map(Broadcasting(Reindex(point_to_x)),cell_to_points)
  cell_to_f = get_array(f)
  cell_to_fxs = lazy_map(evaluate,cell_to_f,cell_to_xs)
  point_to_fxs = lazy_map(Reindex(cell_to_fxs),point_to_cell)
  point_to_fx = lazy_map(getindex,point_to_fxs,point_to_lpoint)
  collect(point_to_fx)          # Collect into a plain array
end

"""
    compute_cell_points_from_vector_of_points(xs::AbstractVector{<:Point},
        trian::Triangulation, domain_style::PhysicalDomain)
"""
function compute_cell_points_from_vector_of_points(xs::AbstractVector{<:Point}, trian::Triangulation, domain_style::PhysicalDomain)
    searchmethod = KDTreeSearch()
    cache1 = _point_to_cell_cache(searchmethod,trian)
    x_to_cell(x) = _point_to_cell!(cache1, x)
    point_to_cell = map(x_to_cell, xs)
    ncells = num_cells(trian)
    cell_to_points, point_to_lpoint = make_inverse_table(point_to_cell, ncells)
    cell_to_xs = lazy_map(Broadcasting(Reindex(xs)), cell_to_points)
    cell_point_xs = CellPoint(cell_to_xs, trian, PhysicalDomain())
end

(a::CellField)(x) = evaluate(a,x)

function evaluate!(cache,f::CellField,x::CellPoint)
  _f, _x = _to_common_domain(f,x)
  cell_field = get_data(_f)
  cell_point = get_data(_x)
  lazy_map(evaluate,cell_field,cell_point)
end

function _to_common_domain(f::CellField,x::CellPoint)
  trian_f = get_triangulation(f)
  trian_x = get_triangulation(x)
  f_on_trian_x = change_domain(f,trian_x,DomainStyle(x))
  f_on_trian_x, x
end

# Gradient

function gradient(a::CellField)
  cell_∇a = lazy_map(Broadcasting(∇),get_data(a))
  if DomainStyle(a) == PhysicalDomain()
    g = cell_∇a
  else
    cell_map = get_cell_map(get_triangulation(a))
    g = lazy_map(Broadcasting(push_∇),cell_∇a,cell_map)
  end
  similar_cell_field(a,g,get_triangulation(a),DomainStyle(a))
end

function DIV(a::CellField)
  DIVa = DIV(get_data(a))
  if DomainStyle(a) == PhysicalDomain()
    @notimplemented
  end
  similar_cell_field(a,DIVa,get_triangulation(a),DomainStyle(a))
end

function ∇∇(a::CellField)
  cell_∇∇a = lazy_map(Broadcasting(∇∇),get_data(a))
  if DomainStyle(a) == PhysicalDomain()
    h = cell_∇∇a
  else
    cell_map = get_cell_map(get_triangulation(a))
    h = lazy_map(Broadcasting(push_∇∇),cell_∇∇a,cell_map)
  end
  similar_cell_field(a,h,get_triangulation(a),DomainStyle(a))
end

# This function has to be removed when ∇⋅∇(a) is implemented
laplacian(a::CellField) = tr(∇∇(a))

# Operations between CellField

function evaluate!(cache,k::Operation,a::CellField...)
  _operate_cellfields(k,a...)
end

function evaluate!(cache,k::Operation,a::Union{Function,CellField}...)
  _operate_cellfields(k,_convert_to_cellfields(a...)...)
end

function evaluate!(cache,k::Operation,a::Union{Number,CellField}...)
  _operate_cellfields(k,_convert_to_cellfields(a...)...)
end

function evaluate!(cache,k::Operation,a::Union{AbstractArray{<:Number},CellField}...)
  _operate_cellfields(k,_convert_to_cellfields(a...)...)
end

# Why julia hangs with this method????
#
#function evaluate!(cache,k::Operation,a::Union{Function,Number,AbstractArray,CellField}...)
#  b = _convert_to_cellfields(a...)
#  _operate_cellfields(k,b...)
#end

struct OperationCellField{DS} <: CellField
  op::Operation
  args::Tuple
  trian::Triangulation
  domain_style::DS
  memo::Dict{Any,Any}
  function OperationCellField(op::Operation,args::CellField...)

    @assert length(args) > 0
    trian = get_triangulation(first(args))
    domain_style = DomainStyle(first(args))
    @check all(i -> DomainStyle(i) == domain_style, args)

    # This is only to catch errors in user code
    # as soon as possible.
    if num_cells(trian) > 0
      @check begin
        pts = _get_cell_points(args...)
        #x = testitem(get_data(pts))
        #f = map(ak -> testitem(get_data(ak)), args)
        #fx = map(fk -> return_value(fk,x), f)
        #r = Fields.BroadcastingFieldOpMap(op.op)(fx...)
        ax = map(i->i(pts),args)
        axi = map(first,ax)
        r = Fields.BroadcastingFieldOpMap(op.op)(axi...)
        true
      end
    end

    new{typeof(domain_style)}(op,args,trian,domain_style,Dict())
  end
end

function _get_cell_points(args::CellField...)
  k = findfirst(i->isa(i,CellState),args)
  if k === nothing
    j = findall(i->isa(i,OperationCellField),args)
    if length(j) == 0
      _get_cell_points(first(args))
    else
      _get_cell_points(args[j]...)
    end
  else
    args[k].points
  end
end

function _get_cell_points(a::CellField)
  trian = get_triangulation(a)
  get_cell_points(trian)
end

function _get_cell_points(a::OperationCellField...)
  b = []
  for ai in a
    for i in ai.args
      push!(b,i)
    end
  end
  _get_cell_points(b...)
end

function _get_cell_points(a::OperationCellField)
  _get_cell_points(a.args...)
end

function get_data(f::OperationCellField)
  a = map(get_data,f.args)
  lazy_map(Broadcasting(f.op),a...)
end
get_triangulation(f::OperationCellField) = f.trian
DomainStyle(::Type{OperationCellField{DS}}) where DS = DS()

function evaluate!(cache,f::OperationCellField,x::CellPoint)
  ax = map(i->i(x),f.args)
  lazy_map(Fields.BroadcastingFieldOpMap(f.op.op),ax...)
end

function change_domain(f::OperationCellField,target_trian::Triangulation,target_domain::DomainStyle)
  args = map(i->change_domain(i,target_trian,target_domain),f.args)
  OperationCellField(f.op,args...)
end

function _operate_cellfields(k::Operation,a...)
  b = _to_common_domain(a...)
  OperationCellField(k,b...)
end

function _convert_to_cellfields(a...)
  a1 = filter(i->isa(i,CellField),a)
  a2 = _to_common_domain(a1...)
  target_domain = DomainStyle(first(a2))
  target_trian = get_triangulation(first(a2))
  map(i->CellField(i,target_trian,target_domain),a)
end

function _to_common_domain(a::CellField...)

  # Find a suitable domain style
  if any( map(i->DomainStyle(i)==ReferenceDomain(),a) )
    target_domain = ReferenceDomain()
  else
    target_domain = PhysicalDomain()
  end

  # Find a suitable triangulation
  msg = """\n
  You are trying to operate CellField objects defined on incompatible triangulations.

  Make sure that all CellField objects are defined on the background triangulation
  or that the number of different sub-triangulations is equal to one.

  For instance:

  - 3 cell fields 2, two them on the same Neumann boundary and the other on the background mesh is OK.

  - 2 cell fields defined on 2 different Neumann boundaries is NOT OK.
  """
  trian_candidates = unique(objectid,map(get_triangulation,a))
  if length(trian_candidates) == 1
    target_trian = first(trian_candidates)
  elseif length(trian_candidates) == 2
    trian_a, trian_b = trian_candidates
    sa_tb = is_change_possible(trian_a,trian_b)
    sb_ta = is_change_possible(trian_b,trian_a)
    if sa_tb && sb_ta
      target_trian = best_target(trian_a,trian_b)
    elseif !sa_tb && sb_ta
      target_trian = trian_a
    elseif sa_tb && !sb_ta
      target_trian = trian_b
    else
      @unreachable msg
    end
  else
    m = """\n
    Cannote operate cellfields defined over more than 2 different
    triangulations at this moment.
    """
    @notimplemented
  end
  map(i->change_domain(i,target_trian,target_domain),a)
end

# Composition (this replaces the @law macro)

Base.:(∘)(f::Function,g::CellField) = Operation(f)(g)
Base.:(∘)(f::Function,g::Tuple{Vararg{CellField}}) = Operation(f)(g...)
Base.:(∘)(f::Function,g::Tuple{Vararg{Union{AbstractArray{<:Number},CellField}}}) = Operation(f)(g...)
Base.:(∘)(f::Function,g::Tuple{Vararg{Union{Function,CellField}}}) = Operation(f)(g...)

# Define some of the well known arithmetic ops

# Unary ops

for op in (:symmetric_part,:inv,:det,:abs,:abs2,:+,:-,:tr,:transpose,:adjoint,:grad2curl,:real,:imag,:conj)
  @eval begin
    ($op)(a::CellField) = Operation($op)(a)
  end
end

# Binary ops

for op in (:inner,:outer,:double_contraction,:+,:-,:*,:cross,:dot,:/)
  @eval begin
    ($op)(a::CellField,b::CellField) = Operation($op)(a,b)
    ($op)(a::CellField,b::Number) = Operation($op)(a,b)
    ($op)(a::Number,b::CellField) = Operation($op)(a,b)
    ($op)(a::CellField,b::Function) = Operation($op)(a,b)
    ($op)(a::Function,b::CellField) = Operation($op)(a,b)
    ($op)(a::CellField,b::AbstractArray{<:Number}) = Operation($op)(a,b)
    ($op)(a::AbstractArray{<:Number},b::CellField) = Operation($op)(a,b)
  end
end

Base.broadcasted(f,a::CellField,b::CellField) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::Number,b::CellField) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::CellField,b::Number) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::Function,b::CellField) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(f,a::CellField,b::Function) = Operation((i,j)->f.(i,j))(a,b)
Base.broadcasted(::typeof(*),::typeof(∇),f::CellField) = Operation(Fields._extract_grad_diag)(∇(f))
Base.broadcasted(::typeof(*),s::Fields.ShiftedNabla,f::CellField) = Operation(Fields._extract_grad_diag)(s(f))

dot(::typeof(∇),f::CellField) = divergence(f)
function (*)(::typeof(∇),f::CellField)
  msg = "Syntax ∇*f has been removed, use ∇⋅f (\\nabla \\cdot f) instead"
  error(msg)
end
outer(::typeof(∇),f::CellField) = gradient(f)
outer(f::CellField,::typeof(∇)) = transpose(gradient(f))
cross(::typeof(∇),f::CellField) = curl(f)

"""
    get_physical_coordinate(trian::Triangulation)

In contrast to get_cell_map, the returned object:
- is a [`CellField`](@ref)
- its gradient is the identity tensor
"""
function get_physical_coordinate(trian::Triangulation)
  CellField(_phys_coord,trian)
end

_phys_coord(x) = x

_phys_coord_grad(x) = one(typeof(outer(x,x)))

gradient(::typeof(_phys_coord)) = _phys_coord_grad

# Skeleton related Operations

function Base.getproperty(x::CellField, sym::Symbol)
  if sym in (:⁺,:plus)
    CellFieldAt{:plus}(x)
  elseif sym in (:⁻, :minus)
    CellFieldAt{:minus}(x)
  else
    getfield(x, sym)
  end
end

function Base.propertynames(x::CellField, private::Bool=false)
  (fieldnames(typeof(x))...,:⁺,:plus,:⁻,:minus)
end

struct CellFieldAt{T,F} <: CellField
  parent::F
  CellFieldAt{T}(parent::CellField) where T = new{T,typeof(parent)}(parent)
end

get_data(f::CellFieldAt) = get_data(f.parent)
get_triangulation(f::CellFieldAt) = get_triangulation(f.parent)
DomainStyle(::Type{CellFieldAt{T,F}}) where {T,F} = DomainStyle(F)
gradient(a::CellFieldAt{P}) where P = CellFieldAt{P}(gradient(a.parent))
∇∇(a::CellFieldAt{P}) where P = CellFieldAt{P}(∇∇(a.parent))
function similar_cell_field(f::CellFieldAt{T},cell_data,trian,ds) where T
  parent = similar_cell_field(f.parent,cell_data,trian,ds)
  CellFieldAt{T}(parent)
end

function CellFieldAt{T}(parent::OperationCellField) where T
  args = map(i->CellFieldAt{T}(i),parent.args)
  OperationCellField(parent.op,args...)
end

function get_normal_vector(trian::Triangulation,cell_normal::SkeletonPair)
  plus = get_normal_vector(trian,cell_normal.plus)
  minus = get_normal_vector(trian,cell_normal.minus)
  SkeletonPair(plus,minus)
end

function get_tangent_vector(trian::Triangulation,cell_tangent::SkeletonPair)
  plus = get_normal_vector(trian,cell_tangent.plus)
  minus = get_normal_vector(trian,cell_tangent.minus)
  SkeletonPair(plus,minus)
end

for op in (:outer,:*,:dot)
  @eval begin
    ($op)(a::CellField,b::SkeletonPair{<:CellField}) = Operation($op)(a,b)
    ($op)(a::SkeletonPair{<:CellField},b::CellField) = Operation($op)(a,b)
  end
end

function evaluate!(cache,k::Operation,a::CellField,b::SkeletonPair{<:CellField})
  plus = k(a.plus,b.plus)
  minus = k(a.minus,b.minus)
  SkeletonPair(plus,minus)
end

function evaluate!(cache,k::Operation,a::SkeletonPair{<:CellField},b::CellField)
  plus = k(a.plus,b.plus)
  minus = k(a.minus,b.minus)
  SkeletonPair(plus,minus)
end

"""
    jump(a::CellField)
    jump(a::SkeletonPair{<:CellField})

Jump operator at interior facets of the supporting `Triangulation`, defined by
`jump`(`a` n) = ⟦`a` n⟧ = `a`⁺n⁺ + `a`⁻n⁻, where n is an oriented normal field
to the interior facets, n⁺ = -n⁻ are the normal pointing into the element on the +
and - side of the facets, and `a`⁺/`a`⁻ are the restrictions of `a` to each
element respectively.
"""
jump(a::CellField) = a.⁺ - a.⁻
jump(a::SkeletonPair{<:CellField}) = a.⁺ + a.⁻ # a.⁻ results from multiplying by n.⁻. Thus we need to sum.

"""
    mean(a::CellField)

Similar to [`jump`](@ref), but for the mean operator `a` ⟶ (`a`⁺ + `a`⁻)/2.
"""
mean(a::CellField) = Operation(_mean)(a.⁺,a.⁻)
_mean(x,y) = 0.5*x + 0.5*y

# This is the fundamental part to make operations on the skeleton work.

for fun in (:change_domain_ref_ref,:change_domain_phys_phys)
  @eval begin

    function $fun(
      a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::SkeletonPair)
      msg = """\n
        It is not possible to use the given CellField on a SkeletonTriangulation.
        Make sure that you are specifying which of the two possible traces,
        either plus (aka ⁺) or minus (aka ⁻) you want to use.
        """
      # If the underlying array is Fill it does not matter if we restrict if from
      # the plus or minus side.
      @check isa(get_array(a.cell_field),Fill)
      plus = $fun(a,ttrian,sglue,tglue.plus)
      plus
    end

    function $fun(
      a::CellFieldAt,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::SkeletonPair)
      if isa(a,CellFieldAt{:plus})
        $fun(a,ttrian,sglue,tglue.plus)
      elseif isa(a,CellFieldAt{:minus})
        $fun(a,ttrian,sglue,tglue.minus)
      else
        @unreachable
      end
    end

    function $fun(
      f::OperationCellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::SkeletonPair)
      args = map(i->$fun(i,ttrian,sglue,tglue),f.args)
      OperationCellField(f.op,args...)
    end

  end
end

# Just to provide more meaningful error messages
function (a::SkeletonPair{<:CellField})(x)
  @unreachable """\n
  You are trying to evaluate a CellField on a mesh skeleton but you have not specified which of the
  two sides i.e. plus (aka ⁺) or minus (aka ⁻) you want to select.

  For instance, if you have extracted the normal vector and the cell points from a SkeletonTriangulation

      x = get_cell_points(strian)
      n = get_normal_vector(strian)

  Evaluating `n(x)` is not allowed. You need to call either `n.⁺(x)` or `n.⁻(x)`.
  """
end

# Interpolable struct
"""
"""
struct KDTreeSearch
  num_nearest_vertices::Int
  function KDTreeSearch(;num_nearest_vertices=1)
    new(num_nearest_vertices)
  end
end

"""
    struct Interpolable{M,A} <: Function
"""
struct Interpolable{M,A} <: Function
  uh::A
  tol::Float64
  searchmethod::M

  """
      Interpolable(uh; tol=1e-6, searchmethod=KDTreeSearch())
  """
  function Interpolable(uh; tol=1e-6, searchmethod=KDTreeSearch())
    new{typeof(searchmethod),typeof(uh)}(uh, tol,searchmethod)
  end
end

(a::Interpolable)(x) = evaluate(a,x)
evaluate!(cache,a::Interpolable,x::Point) = evaluate!(cache,a.uh,x)
return_cache(f::CellField,x::Point) = return_cache(Interpolable(f),x)

function return_cache(a::Interpolable,x::Point)
  f = a.uh
  trian = get_triangulation(f)
  cache1 = _point_to_cell_cache(a.searchmethod,trian)

  cell_f = get_array(f)
  cell_f_cache = array_cache(cell_f)
  cf = testitem(cell_f)
  f_cache = return_cache(cf,x)
  cache2 = cell_f_cache, f_cache, cell_f, f

  return cache1,cache2
end

function _point_to_cell_cache(searchmethod::KDTreeSearch,trian::Triangulation)
  model = get_active_model(trian)
  topo = get_grid_topology(model)
  vertex_coordinates = Geometry.get_vertex_coordinates(topo)
  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
  D = num_cell_dims(trian)
  vertex_to_cells = get_faces(topo, 0, D)
  cell_to_ctype = get_cell_type(trian)
  ctype_to_reffe = get_reffes(trian)
  ctype_to_polytope = map(get_polytope, ctype_to_reffe)
  cell_map = get_cell_map(trian)
  table_cache = array_cache(vertex_to_cells)
  cache1 = searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache
end
