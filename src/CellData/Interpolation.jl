
# This file implements code to evaluate CellFields on arbitrary points, based 
# on tree searches.

# Interpolable struct
struct KDTreeSearch{T}
  num_nearest_vertices::Int
  tol::T
  function KDTreeSearch(; num_nearest_vertices=1, tol=1.e-10)
    T = typeof(tol)
    new{T}(num_nearest_vertices, tol)
  end
end

struct Interpolable{M,A} <: Function
  uh::A
  tol::Float64
  searchmethod::M
  function Interpolable(uh; tol=1e-10, searchmethod=KDTreeSearch(; tol=tol))
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

  return cache1, cache2
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

function evaluate!(cache,f::CellField,point_to_x::AbstractVector{<:Point})
  cache1, cache2 = cache
  searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map = cache1
  cell_f_cache, f_cache, cell_f, f₀ = cache2
  @check f === f₀ "Wrong cache"

  ncells = length(cell_map)
  x_to_cell(x) = _point_to_cell!(cache1,x)
  point_to_cell = map(x_to_cell,point_to_x)
  cell_to_points, point_to_lpoint = make_inverse_table(point_to_cell,ncells)
  cell_to_xs = lazy_map(Broadcasting(Reindex(point_to_x)),cell_to_points)
  cell_to_f = get_array(f)
  cell_to_fxs = lazy_map(evaluate,cell_to_f,cell_to_xs)
  point_to_fxs = lazy_map(Reindex(cell_to_fxs),point_to_cell)
  point_to_fx = lazy_map(getindex,point_to_fxs,point_to_lpoint)
  return collect(point_to_fx) # Collect into a plain array
end

function compute_cell_points_from_vector_of_points(
  xs::AbstractVector{<:Point}, trian::Triangulation, domain_style::PhysicalDomain
)
  searchmethod = KDTreeSearch()
  cache = _point_to_cell_cache(searchmethod,trian)
  x_to_cell(x) = _point_to_cell!(cache, x)
  point_to_cell = map(x_to_cell, xs)
  ncells = num_cells(trian)
  cell_to_points, point_to_lpoint = make_inverse_table(point_to_cell, ncells)
  cell_to_xs = lazy_map(Broadcasting(Reindex(xs)), cell_to_points)
  return CellPoint(cell_to_xs, trian, PhysicalDomain())
end

# Helpers

function _point_to_cell_cache(searchmethod::KDTreeSearch,trian::Triangulation)
  model = get_active_model(trian)
  topo = get_grid_topology(model)
  if num_nodes(model) == num_vertices(model)
    # Non-periodic case
    vertex_coordinates = Geometry.get_vertex_coordinates(topo)
    vertex_to_cells = get_faces(topo, 0, num_cell_dims(trian))
  else
    # Periodic case
    vertex_coordinates = collect1d(Geometry.get_node_coordinates(model))
    cell_to_vertices = Table(get_cell_node_ids(model))
    vertex_to_cells = Arrays.inverse_table(cell_to_vertices, num_nodes(model))
  end
  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
  cell_to_ctype = get_cell_type(trian)
  ctype_to_polytope = get_polytopes(trian)
  cell_map = get_cell_map(trian)
  table_cache = array_cache(vertex_to_cells)
  return searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache
end

function _point_to_cell!(cache, x::Point)
  searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache

  function cell_distance(cell::Integer)
    ctype = cell_to_ctype[cell]
    polytope = ctype_to_polytope[ctype]
    cmap = cell_map[cell]
    inv_cmap = inverse_map(cmap)
    return distance(polytope, inv_cmap, x)
  end

  # Find the nearest vertices to the point `x` in the triangulation
  vertices, distances = knn(kdtree, get_array(ForwardDiff.value(x)), searchmethod.num_nearest_vertices, true)

  T = eltype(distances)
  tol = max(1000*eps(T), T(searchmethod.tol))
  for vertex in vertices

    # Find all neighbouring cells
    cells = getindex!(table_cache,vertex_to_cells,vertex)
    @assert !isempty(cells)

    # Calculate the distance from the point to all the neighbor cells. Without
    # round-off, and with non-overlapping cells, the distance would be
    # negative for exactly one cell and positive for all other ones. Due
    # to round-off, the distance can be slightly negative or slightly
    # positive for points near cell boundaries, in particular near
    # vertices. In this case, choose the cell with the smallest
    # distance, and check that the distance (if positive) is at most at
    # round-off level.

    cell = zero(eltype(cells))
    dist = T(Inf)
    for jcell in cells
      jdist = cell_distance(jcell)
      if jdist < dist
        cell = jcell
        dist = jdist
      end
    end

    (dist < tol) && return cell
  end

  # Output error message if cell not found
  @check false "Point $x is not inside any active cell"
end

"""
   dist = distance(polytope::ExtrusionPolytope,inv_cmap::Field,x::Point)

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
