
# Descriptor of a cartesian grid

"""
    struct CartesianDescriptor{D,T,F<:Function}
      origin::Point{D,T}
      sizes::NTuple{D,T}
      partition::NTuple{D,Int}
      map::F
    end

Struct that stores the data defining a Cartesian grid.
"""
struct CartesianDescriptor{D,T,F<:Function} <: GridapType
  origin::Point{D,T}
  sizes::NTuple{D,T}
  partition::NTuple{D,Int}
  map::F
  isperiodic::NTuple{D,Bool}
  @doc """
      CartesianDescriptor(
        origin::Point{D},
        sizes::NTuple{D},
        partition;
        map::Function=identity,
        isperiodic::NTuple{D,Bool}=tfill(false,Val{D})) where D

  `partition` is a 1D indexable collection of arbitrary type.
  """
  function CartesianDescriptor(
    origin::Point{D},
    sizes::NTuple{D},
    partition;
    map::Function=identity,
    isperiodic::NTuple{D,Bool}=tfill(false,Val{D}())) where D

    for i in 1:D
      if isperiodic[i]
        @assert (partition[i] > 2) "A minimum of 3 elements is required in any
          periodic direction"
      end
    end

    T = eltype(sizes)
    F = typeof(map)
    new{D,T,F}(origin,sizes,Tuple(partition),map,isperiodic)
  end
end

"""
    CartesianDescriptor(
      domain,
      partition;
      map::Function=identity,
      isperiodic::NTuple{D,Bool}=tfill(false,Val{D}))

`domain` and `partition` are 1D indexable collections of arbitrary type.
"""
function CartesianDescriptor(
  domain,
  partition;
  map::Function=identity,
  isperiodic::NTuple=tfill(false,Val{length(partition)}()))

  D = length(partition)
  limits = [(domain[2*d-1],domain[2*d]) for d in 1:D]
  sizes = Tuple([(limits[d][2]-limits[d][1])/partition[d] for d in 1:D])
  origin = Point([ limits[d][1] for d in 1:D]...)
  CartesianDescriptor(origin,sizes,partition;map=map,isperiodic=isperiodic)
end

"""
    CartesianDescriptor(
      pmin::Point{D},
      pmax::Point{D},
      partition;
      map::Function=identity,
      isperiodic::NTuple{D,Bool}=tfill(false,Val{D})) where D

`partition` is a 1D indexable collection of arbitrary type.
"""
function CartesianDescriptor(
  pmin::Point{D},
  pmax::Point{D},
  partition;
  map::Function=identity,
  isperiodic::NTuple{D,Bool}=tfill(false,Val{D}())) where D

  T = eltype(pmin)
  domain = zeros(T,2*D)
  for d in 1:D
    domain[2*(d-1)+1] = pmin[d]
    domain[2*(d-1)+2] = pmax[d]
  end
  CartesianDescriptor(domain,partition;map=map,isperiodic=isperiodic)
end

# Deprecated signatures for backwards compatibility

@deprecate CartesianDescriptor(origin::Point{D}, sizes::NTuple{D}, partition, map::Function) where D CartesianDescriptor(origin, sizes, partition; map=map)

@deprecate CartesianDescriptor(domain,partition,map::Function) CartesianDescriptor(domain,partition;map=map)

@deprecate CartesianDescriptor(pmin::Point{D}, pmax::Point{D}, partition, map::Function) where D CartesianDescriptor(pmin, pmax, partition; map=map)

# Coordinates

struct CartesianCoordinates{D,T,F} <: AbstractArray{Point{D,T},D}
  data::CartesianDescriptor{D,T,F}
end

Base.size(a::CartesianCoordinates) = a.data.partition .+ 1

Base.IndexStyle(::Type{<:CartesianCoordinates}) = IndexCartesian()

function Base.getindex(a::CartesianCoordinates{D,T}, I::Vararg{Int,D}) where {D,T}
  p = zero(mutable(Point{D,T}))
  x0 = a.data.origin
  dx = a.data.sizes
  @inbounds for d in 1:D
    p[d] =  x0[d] + (I[d]-1)*dx[d]
  end
  a.data.map(Point(p))
end

# Cell nodes

struct CartesianCellNodes{D} <: AbstractArray{Vector{Int},D}
  partition::NTuple{D,Int}
  function CartesianCellNodes(partition)
    D = length(partition)
    new{D}(partition)
  end
end

Base.size(a::CartesianCellNodes) = a.partition

Base.IndexStyle(::Type{<:CartesianCellNodes}) = IndexCartesian()

function array_cache(a::CartesianCellNodes{D}) where D
  zeros(Int,2^D)
end

@inline function getindex!(cache,a::CartesianCellNodes,i::Integer)
  cis = CartesianIndices(size(a))
  ci = cis[i]
  getindex!(cache,a,Tuple(ci)...)
end

@inline function getindex!(v,a::CartesianCellNodes{1},i::Integer)
  v[1] = i
  v[2] = i+1
  v
end

@inline function getindex!(v,a::CartesianCellNodes{D},i::Integer...) where D
  nodes = LinearIndices(size(a).+1)
  lnodes = CartesianIndices(tfill(2,Val{D}()))
  j = i .- 1
  for (i,lnode) in enumerate(lnodes)
    k = Tuple(lnode) .+ j
    v[i] = nodes[ k... ]
  end
  v
end

function Base.getindex(a::CartesianCellNodes,i::Integer...)
  cache = array_cache(a)
  getindex!(cache,a,i...)
end

# Cartesian Grid

"""
    struct CartesianGrid{D,T,F} <: Grid{D,D}
      # private fields
    end
"""
struct CartesianGrid{D,T,F} <: Grid{D,D}
  node_coords::CartesianCoordinates{D,T,F}
  cell_nodes::CartesianCellNodes{D}
  cell_type::Fill{Int8,1,Tuple{Base.OneTo{Int}}}
  memo::Dict
  @doc """
      CartesianGrid(desc::CartesianDescriptor)
  """
  function CartesianGrid(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
    node_coords = CartesianCoordinates(desc)
    cell_nodes = CartesianCellNodes(desc.partition)
    cell_type = Fill(Int8(1),length(cell_nodes))
    new{D,T,F}(node_coords,cell_nodes,cell_type,Dict())
  end
end

get_memo(a::CartesianGrid) = a.memo

OrientationStyle(::Type{<:CartesianGrid}) = Val{true}()

"""
    get_cartesian_descriptor(grid::CartesianGrid)

Get the descriptor of the Cartesian grid
"""
function get_cartesian_descriptor(a::CartesianGrid)
  a.node_coords.data
end

get_node_coordinates(g::CartesianGrid) = g.node_coords

get_cell_type(g::CartesianGrid) = g.cell_type

get_cell_nodes(g::CartesianGrid) = g.cell_nodes

function get_reffes(g::CartesianGrid{D}) where D
  p = Polytope(tfill(HEX_AXIS,Val{D}()))
  order = 1
  reffe = LagrangianRefFE(Float64,p,order)
  [reffe,]
end

get_reffes(g::CartesianGrid{1}) = [SEG2,]

get_reffes(g::CartesianGrid{2}) = [QUAD4,]

get_reffes(g::CartesianGrid{3}) = [HEX8,]

"""
    CartesianGrid(args...;kwargs...)

Same args needed to construct a `CartesianDescriptor`
"""
function CartesianGrid(args...;kwargs...)
  desc = CartesianDescriptor(args...;kwargs...)
  CartesianGrid(desc)
end

# Cell map

struct CartesianMap{D,T,L} <: AbstractArray{AffineMap{D,T,L},D}
  data::CartesianDescriptor{D,T,typeof(identity)}
  function CartesianMap(des::CartesianDescriptor{D,T}) where {D,T}
    L = D*D
    new{D,T,L}(des)
  end
end

Base.size(a::CartesianMap) = a.data.partition

Base.IndexStyle(::Type{<:CartesianMap}) = IndexCartesian()

function Base.getindex(a::CartesianMap{D,T},I::Vararg{Int,D}) where {D,T}
  p = zero(mutable(Point{D,T}))
  x0 = a.data.origin
  dx = a.data.sizes
  @inbounds for d in 1:D
    p[d] =  x0[d] + (I[d]-1)*dx[d]
  end
  origin =  Point(p)
  jacobian = diagonal_tensor(VectorValue(dx))
  AffineMap(jacobian,origin)
end

function field_array_gradient(a::CartesianMap)
  dx = a.data.sizes
  jacobian = diagonal_tensor(VectorValue(dx))
  j = AffineMapGrad(jacobian)
  Fill(j,length(a))
end

function field_array_gradient(a::Reindexed{T,N,A}) where {T,N,A<:CartesianMap}
  g = field_array_gradient(a.i_to_v)
  reindex(g,a.j_to_i)
end

function compute_cell_map(grid::CartesianGrid{D,T,typeof(identity)} where {D,T})
  array = CartesianMap(grid.node_coords.data)
  GenericCellMap(array)
end

# Cartesian grid topology with periodic BC

function _cartesian_grid_topology_with_periodic_bcs(grid::UnstructuredGrid,
  isperiodic::NTuple,
  partition)

  cell_to_vertices, vertex_to_node =
    _generate_cell_to_vertices_from_grid(grid, isperiodic, partition)
  _generate_grid_topology_from_grid(grid,cell_to_vertices,vertex_to_node)
end

function _generate_cell_to_vertices_from_grid(grid::UnstructuredGrid,
  isperiodic::NTuple, partition)

  if is_first_order(grid)
    nodes = get_cell_nodes(grid)
    cell_to_vertices = copy(nodes)

    nnodes = num_nodes(grid)
    num_nodes_x_dir = [partition[i]+1 for i in 1:length(partition)]
    point_to_isperiodic, slave_point_to_point, slave_point_to_master_point =
      _generate_slave_to_master_point(num_nodes_x_dir,isperiodic, nnodes)

    vertex_to_point = findall( .! point_to_isperiodic)
    point_to_vertex = fill(-1,length(point_to_isperiodic))
    point_to_vertex[vertex_to_point] = 1:length(vertex_to_point)
    point_to_vertex[slave_point_to_point] = point_to_vertex[slave_point_to_master_point]

    cell_to_vertices = Table(LocalToGlobalArray(nodes,point_to_vertex))

    vertex_to_node = vertex_to_point
    node_to_vertex = point_to_vertex
  else
    @notimplemented
  end
  (cell_to_vertices,vertex_to_node, node_to_vertex)
end

function _generate_slave_to_master_point(num_nodes_x_dir::Vector{Int},
  isperiodic::NTuple, num_nodes::Int)

  point_to_isperiodic = fill(false,num_nodes)

  slave_ijk_bounds = Array{Any,1}(undef,length(isperiodic))
  master_ijk_bounds = Array{Any,1}(undef,length(isperiodic))

  linear_indices = LinearIndices(Tuple(num_nodes_x_dir))
  periodic_dirs = findall(x->x==true, isperiodic)
  for periodic_dir in periodic_dirs
    for i in 1:length(isperiodic)
      if i == periodic_dir
        slave_ijk_bounds[i] = num_nodes_x_dir[i]
      else
        slave_ijk_bounds[i] = 1:num_nodes_x_dir[i]
      end
    end
    slave_ijk_set = CartesianIndices(Tuple(slave_ijk_bounds))
    point_to_isperiodic[linear_indices[slave_ijk_set]] .= true
  end

  slave_point_to_point = findall( point_to_isperiodic)
  slave_point_to_master_point = Array{Int,1}(undef,length(slave_point_to_point))

  cartesian_indices = CartesianIndices(Tuple(num_nodes_x_dir))
  for (i,point) in enumerate(slave_point_to_point)
    ijk = collect(cartesian_indices[point].I)
    for i in periodic_dirs
      if ijk[i] == num_nodes_x_dir[i]
        ijk[i] = 1
      end
    end

    master_point_ijk = CartesianIndex(Tuple(ijk))
    slave_point_to_master_point[i] = linear_indices[master_point_ijk]
  end

  point_to_isperiodic, slave_point_to_point, slave_point_to_master_point
end
