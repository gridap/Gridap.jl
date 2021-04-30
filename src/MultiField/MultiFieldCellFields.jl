
struct MultiFieldCellField{DS<:DomainStyle} <: CellField
  single_fields::Vector{<:CellField}
  trian::Triangulation
  domain_style::DS

  function MultiFieldCellField(single_fields::Vector{<:CellField})
    @assert length(single_fields) > 0
    f1 = first(single_fields)

    # Find a suitable domain style
    if any( map(i->DomainStyle(i)==ReferenceDomain(),single_fields) )
      domain_style = ReferenceDomain()
    else
      domain_style = PhysicalDomain()
    end

    trian = get_triangulation(f1)
    #@check all(map(i-> get_triangulation(i) === trian,single_fields))

    new{typeof(domain_style)}(single_fields,trian,domain_style)
  end
end

function CellData.get_data(f::MultiFieldCellField)
  s = """
  Function get_data is not implemented for MultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separatelly.

  If ever implement this, evaluating a `MultiFieldCellField` directly would provide,
  at each evaluation point, a tuple with the value of the different fields.
  """
  @notimplemented s
end

CellData.get_triangulation(f::MultiFieldCellField) = f.trian
CellData.DomainStyle(::Type{MultiFieldCellField{DS}}) where DS = DS()
num_fields(a::MultiFieldCellField) = length(a.single_fields)
Base.getindex(a::MultiFieldCellField,i::Integer) = a.single_fields[i]
Base.iterate(a::MultiFieldCellField)  = iterate(a.single_fields)
Base.iterate(a::MultiFieldCellField,state)  = iterate(a.single_fields,state)

# Evaluation of CellFields

# This code lives here (instead of in the CellData module) because we
# need access to MultiField functionality.

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

function return_cache(f::CellField,x::Point)
  trian = get_triangulation(f)
  topo = GridTopology(trian)    # Note: this is expensive
  vertex_coordinates = Geometry.get_vertex_coordinates(topo)
  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
  D = num_cell_dims(trian)
  vertex_to_cells = get_faces(topo,0,D)
  cell_to_ctype = get_cell_type(trian)
  ctype_to_reffe = get_reffes(trian)
  ctype_to_polytope = map(get_polytope, ctype_to_reffe)
  cell_map = get_cell_map(trian)
  cache1 = kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map

  cell_f = get_array(f)
  cell_f_cache = array_cache(cell_f)
  cf = testitem(cell_f)
  f_cache = return_cache(cf,x)
  cache2 = cell_f_cache, f_cache, cell_f, f

  return cache1,cache2
end

function point_to_cell!(cache,f::CellField,x::Point)
  kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map = cache

  # Find nearest vertex
  id,dist = nn(kdtree, SVector(Tuple(x)))

  # Find all neighbouring cells
  cells = vertex_to_cells[id]
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
  # Ensure the point is inside one of the cells, up to round-off errors
  @check dist ≤ 1000eps(T) "Point is not inside any cell"

  return cell
end

function evaluate!(cache,f::CellField,x::Point)
  cache1,cache2 = cache
  cell_f_cache, f_cache, cell_f, f₀ = cache2
  @check f === f₀ "Wrong cache"

  cell = point_to_cell!(cache1,f,x)
  cf = getindex!(cell_f_cache, cell_f, cell)
  fx = evaluate!(f_cache, cf, x)
  return fx
end

return_cache(f::CellField,xs::AbstractVector{<:Point}) = return_cache(f,testitem(xs))

# # Simple version:
# function evaluate!(cache,f::CellField,xs::AbstractVector{<:Point})
#   return map(x->evaluate!(cache,f,x), xs)
# end

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
  cache1,cache2 = return_cache(f,testitem(point_to_x))
  kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map = cache1
  cell_f_cache, f_cache, cell_f, f₀ = cache2
  @check f === f₀ "Wrong cache"

  ncells = length(cell_map)
  x_to_cell(x) = point_to_cell!(cache1,f,x)
  point_to_cell = map(x_to_cell,point_to_x)
  cell_to_points,point_to_lpoint = make_inverse_table(point_to_cell,ncells)
  cell_to_xs = lazy_map(Broadcasting(Reindex(point_to_x)),cell_to_points)
  cell_to_f = get_array(f)
  cell_to_fxs = lazy_map(evaluate,cell_to_f,cell_to_xs)
  point_to_fxs = lazy_map(Reindex(cell_to_fxs),point_to_cell)
  point_to_fx = lazy_map(getindex,point_to_fxs,point_to_lpoint)
  collect(point_to_fx)          # Collect into a plain array
end
