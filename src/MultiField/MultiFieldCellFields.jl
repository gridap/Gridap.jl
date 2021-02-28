
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

mutable struct CellFieldCache
  kdtree::Union{Nothing, KDTree{SVector{D,T},Euclidean,T} where {D,T}}
end

# TODO: Should this return the cache already in `f` instead of
# creating a new one?
return_cache(f::CellField,x::Point) = CellFieldCache(nothing)

function evaluate!(cache,f::CellField,x::Point)
  cache::CellFieldCache

  # Examine triangulation
  trian = get_triangulation(f)
  node_coordinates = get_node_coordinates(trian)
  @assert !isempty(node_coordinates)
  sample_coord = node_coordinates[1]
  # TODO: call Gridap.ReferenceFEs.num_cell_dims instead?
  D = length(sample_coord)
  T = eltype(sample_coord)
  @assert D > 0

  # Find polytope
  reffes = trian.reffes
  @assert length(reffes) == 1
  reffe = reffes[1].reffe
  polytope = reffe.polytope
  extrusion = polytope.extrusion
  @assert length(extrusion) == D

  if all(e == TET_AXIS for e in extrusion)

    # Get k-d tree for mesh from cache
    if cache.kdtree === nothing
      cache.kdtree = KDTree(map(nc -> SVector(Tuple(nc)), node_coordinates))
    end
    kdtree = cache.kdtree::KDTree{SVector{D,T},Euclidean,T}

    # Find all neighbouring cells
    id,dist = nn(kdtree, SVector(Tuple(x)))
    cells = Int[]
    # TODO: Use better algorithm
    #
    # Francesc Verdugo @fverdugo Dec 13 2020 12:46
    # vertex_to_cells = get_faces(get_grid_topology(model),0,num_cell_dims(model))
    # This will provide cells around each vertex (vertex == node for
    # linear meshes)
    #
    # Francesc Verdugo @fverdugo 14:07
    # You have to start with a DiscreteModel, then
    # topo=get_grid_topology(model)
    # vertex_to_cells = get_faces(topo,0,num_cell_dims(model))
    for (c,ids) in enumerate(trian.cell_node_ids)
      id ∈ ids && push!(cells, c)
    end
    @assert !isempty(cells)

    # Calculate barycentric coordinates for a point x in the given cell
    function bary(cell, x)
      local ids = trian.cell_node_ids[cell]
      @assert length(ids) == D+1 "simplex has correct number of vertices"
      local ys = zero(SVector{D+1,SVector{D,T}})
      for (n,id) in enumerate(ids)
        local y = SVector(Tuple(node_coordinates[id]))
        ys = setindex(ys, y, n)
      end
      local λ = cartesian2barycentric(ys, SVector{D,T}(Tuple(x)))
      return λ
    end

    # Calculate the distance from the point to all the cells. Without
    # round-off, and with non-overlapping cells, the distance would be
    # negative for exactly one cell and positive for all other ones.
    # Due to round-off, the distance can be slightly negative or
    # slightly positive for points near cell boundaries, in particular
    # near vertices. In this case, choose the cell with the smallest
    # distance, and check that the distance (if positive) is at most
    # at round-off level.
    mindist = T(Inf)
    minc,minλ = 0, zero(SVector{D+1,T})
    for c in cells
      λ = bary(c, x)
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

    # TODO: Actually evaluate basis functions
    # Francesc Verdugo @fverdugo Dec 09 2020 13:06
    # for 2., take a look at the functions in this file https://github.com/gridap/Gridap.jl/blob/master/src/Geometry/FaceLabelings.jl
    # In particular this one:
    # https://github.com/gridap/Gridap.jl/blob/4342137dfbe635867f39e64ecb01e8dbbc5fa6e7/src/Geometry/FaceLabelings.jl#L289
    values = get_cell_dof_values(f)[cell]
    fx = sum(λ[n] * values[n] for n in 1:D+1)

    # you can just evaluate the (finite element) solution at this
    # point at this cell.

    # TODO: If you define the finite element space in the physical
    # space, i.e., DomainStyle being PhysicalDomain(), that is all.
    # However, if DomainStyle equal to ReferenceDomain() one needs the
    # inverse geometrical map that takes the physical cell and maps it
    # to the reference one (in which we evaluate shape functions)

    # TODO: are these lines useful?
    # v = get_cell_shapefuns(V)
    # u = get_cell_shapefuns_trial(U)
    # cellmat = get_array( ∫( ∇(u)⋅∇(v) )dΩ )

  elseif all(e == HEX_AXIS for e in extrusion)

    # TOOD
    @assert false "n-cube not implemented"

  else

    @assert false "Unsupported polytope"

  end

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
