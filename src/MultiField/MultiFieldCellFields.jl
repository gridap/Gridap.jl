
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
    @check all(map(i-> get_triangulation(i) === trian,single_fields))

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

function return_cache(f::CellField,x::Point)
  trian = get_triangulation(f)
  node_coordinates = get_node_coordinates(trian)
  node_coordinates = reshape(node_coordinates, length(node_coordinates))
  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), node_coordinates))
  # TODO: This does not work for SingleFieldFEFunction
  cell_f = get_array(f)
  c1 = array_cache(cell_f)
  f = testitem(cell_f)
  c2 = return_cache(f,x)
  return kdtree, c1, c2, cell_f, f
end

function evaluate!(cache,f::CellField,x::Point)
  kdtree, c1, c2, cell_g, g = cache
  if f === g
    cell_f = cell_g
  else
    cell_f = get_array(f)
  end

  # Find nearest vertex
  id,dist = nn(kdtree, SVector(Tuple(x)))
  @show id
  T = typeof(dist)

  # Find all neighbouring cells
  trian = get_triangulation(f)
  @show typeof(trian)
  node_coordinates = get_node_coordinates(trian)
  node_coordinates = reshape(node_coordinates, length(node_coordinates))
  model = trian   # this is wrong
  # topo = get_grid_topology(model)
  topo = GridTopology(model)
  D = num_cell_dims(trian)
  vertex_to_cells = get_faces(topo,0,D)
  cells = vertex_to_cells[id]
  @show cells

  # Find polytope
  # TODO: this only works for simplices, not for hypercubes
  # type CartesianGrid has no field reffes
  @show typeof(trian)
  reffes = trian.reffes
  @assert length(reffes) == 1
  reffe = reffes[1].reffe
  polytope = reffe.polytope
  extrusion = polytope.extrusion
  @assert length(extrusion) == D

  if all(e == TET_AXIS for e in extrusion)

    # Calculate the distance from the point to all the cells. Without
    # round-off, and with non-overlapping cells, the distance would be
    # positive for exactly one cell and positive for all other ones.
    # Due to round-off, the distance can be slightly negative or
    # slightly positive for points near cell boundaries, in particular
    # near vertices. In this case, choose the cell with the largest
    # distance, and check that the distance (if negative) is at most
    # at round-off level.
    maxdist = -T(Inf)
    maxcell = 0
    for cell in cells
      @show "for" cell
      # Calculate barycentric coordinates λ for x in cell c
      ids = trian.cell_node_ids[cell]
      @show "for" ids
      N = length(ids)
      # @show N D+1
      @assert N == D+1 "simplex has correct number of vertices"
      ys = zero(SVector{N,SVector{D,T}})
      for (n,id) in enumerate(ids)
        y = SVector(Tuple(node_coordinates[id]))
        ys = setindex(ys, y, n)
      end
      @show "for" ys
      λ = cartesian2barycentric(ys, SVector{D,T}(Tuple(x)))
      @show "for" λ
      # dist is the distance to the closest face. If positive, the point
      # is inside this cell.
      dist = minimum(λ)
      @show "for" dist
      if dist > maxdist
        maxcell = cell
        maxdist = dist
      end
    end
    @assert maxcell > 0
    cell = maxcell
    dist = maxdist
  
    # Ensure the point is inside one of the cells, up to round-off errors
    # TODO: Add proper run-time check
    @show cell dist
    @assert dist ≥ -1000 * eps(T) "Point is not inside a cell"

  elseif all(e == HEX_AXIS for e in extrusion)

    maxdist = -T(Inf)
    maxcell = 0
    for cell in cells
      ids = trian.cell_node_ids[cell]
      N = length(ids)
      @assert N == 2^D "hypercube has correct number of vertices"
      dist = T(Inf)
      for d in 1:D
        for f in 1:2
          face = [i+1 for i in 0:2^D-1 if (i & (1<<(D-d)) == 0) == (f==1)]
          yc = ids[face[1]][d]
          @assert all(ids[j][d] == y for j in face)
          dist = min(dist, f == 1 ? x[d] - y : y - x[d])
        end
      end
      if dist > maxdist
        maxcell = cell
        maxdist = dist
      end
    end
    @assert maxcell > 0
    cell = maxcell
    dist = maxdist
  
    # Ensure the point is inside one of the cells, up to round-off errors
    # TODO: Add proper run-time check
    @show cell dist
    @assert dist ≥ -1000 * eps(T) "Point is not inside a cell"

  else

    @notimplemented "Need either tet or hex mesh"

  end

  f = getindex!(c1,cell_f,cell)
  fx = evaluate!(c2,f,x)
  return fx
end

"""
    λ = cartesian2barycentric(s, p)

Convert Cartesian to barycentric coordinates.

# Arguments
- `s`: Simplex vertices in Cartesian coordinates. `s` has `N ≤ D + 1`
  vertices in `D` dimensions.
- `p`: Point in Cartesian coordinates

# Result
- `λ`: Point in barycentric coordinates
"""
function cartesian2barycentric(s::SMatrix{N,D,T}, p::SVector{D,T}) where {N,D,T}
  @assert N ≤ D + 1
  # Algorithm as described on
  # <https://en.wikipedia.org/wiki/Barycentric_coordinate_system>,
  # section "Conversion between barycentric and Cartesian coordinates"
  A = SMatrix{D+1,N}(i == D+1 ? T(1) : s[j, i] for i in 1:D+1, j in 1:N)
  b = SVector{D+1}(p..., T(1))
  return A \ b
end

function cartesian2barycentric(s::SVector{N,SVector{D,T}},
                               p::SVector{D,T}) where {D,N,T}
  return cartesian2barycentric(SMatrix{N,D,T}(s[i][j] for i in 1:N, j in 1:D),
                               p)
end

function barycentric2cartesian(s::SVector{N,SVector{D,T}},
                               λ::SVector{N,T}) where {D,N,T}
  @assert N ≤ D + 1
  return SVector{D,T}(sum(s[i][j] * λ[i] for i in 1:N) for j in 1:D)
end
