
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
  cell_f = get_array(f)
  c1 = array_cache(cell_f)
  f = testitem(cell_f)
  @show typeof(f) typeof(x) typeof(cell_f)
  c2 = return_cache(f,x)
  return kdtree, c1, c2, cell_f, f
end

# function evaluate!(cache,f::CellField,x::Point)
function evaluate!(cache,f::SingleFieldFEFunction,x::Point)
  kdtree, c1, c2, cell_g, g = cache
  if f === g
    cell_f = cell_g
  else
    cell_f = get_array(f)
  end

  # Find nearest vertex
  id,dist = nn(kdtree, SVector(Tuple(x)))

  # Find all neighbouring cells
  trian = get_triangulation(f)
  # topo = get_grid_topology(trian)
  topo = GridTopology(trian)
  vertex_to_cells = get_faces(topo,0,num_cell_dims(model))
  cells = vertex_to_cells[id]

  # Calculate the distance from the point to all the cells. Without
  # round-off, and with non-overlapping cells, the distance would be
  # negative for exactly one cell and positive for all other ones. Due
  # to round-off, the distance can be slightly negative or slightly
  # positive for points near cell boundaries, in particular near
  # vertices. In this case, choose the cell with the smallest
  # distance, and check that the distance (if positive) is at most at
  # round-off level.
  mindist = T(Inf)
  minc,minλ = 0,zero(SVector{D+1,T})
  for c in cells
    # Calculate barycentric coordinates λ for x in cell c
    ids = trian.cell_node_ids[c]
    @assert length(ids) == D+1 "simplex has correct number of vertices"
    ys = zero(SVector{D+1,SVector{D,T}})
    for (n,id) in enumerate(ids)
      y = SVector(Tuple(node_coordinates[id]))
      ys = setindex(ys, y, n)
    end
    λ = cartesian2barycentric(ys, SVector{D,T}(Tuple(x)))
    dist = max(-minimum(λ), maximum(λ) - 1)
    if dist < mindist
      minc,minλ = c,λ
      mindist = dist
    end
  end
  @assert 1 ≤ minc
  cell,λ = minc,minλ
  dist = mindist

  # Ensure the point is inside one of the cells, up to round-off errors
  # TODO: Add proper run-time check
  @assert dist ≤ 1000 * eps(T) "Point is not inside a cell"

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
