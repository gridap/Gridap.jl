
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

function Base.getindex(a::CartesianCoordinates{D,T}, I::Vararg{Integer,D}) where {D,T}
  p = zero(Mutable(Point{D,T}))
  x0 = a.data.origin
  dx = a.data.sizes
  @inbounds for d in 1:D
    p[d] =  x0[d] + (I[d]-1)*dx[d]
  end
  a.data.map(Point(p))
end

# Cell nodes

struct CartesianCellNodes{D} <: AbstractArray{Vector{Int32},D}
  partition::NTuple{D,Int}
  function CartesianCellNodes(partition)
    D = length(partition)
    new{D}(partition)
  end
end

Base.size(a::CartesianCellNodes) = a.partition

Base.IndexStyle(::Type{<:CartesianCellNodes}) = IndexCartesian()

function array_cache(a::CartesianCellNodes{D}) where D
  zeros(Int32,2^D)
end

#@inline function getindex!(cache,a::CartesianCellNodes,i::Integer)
#  cis = CartesianIndices(size(a))
#  ci = cis[i]
#  getindex!(cache,a,Tuple(ci)...)
#end
#
#@inline function getindex!(v,a::CartesianCellNodes{1},i::Integer)
#  v[1] = i
#  v[2] = i+1
#  v
#end

@inline function getindex!(v,a::CartesianCellNodes{D},i::Vararg{Integer,D}) where D
  nodes = LinearIndices(size(a).+1)
  lnodes = CartesianIndices(tfill(Int32(2),Val{D}()))
  j = i .- 1
  for (i,lnode) in enumerate(lnodes)
    k = Tuple(lnode) .+ j
    v[i] = nodes[ k... ]
  end
  v
end

function Base.getindex(a::CartesianCellNodes{D},i::Vararg{Integer,D}) where D
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
  @doc """
      CartesianGrid(desc::CartesianDescriptor)
  """
  function CartesianGrid(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
    node_coords = CartesianCoordinates(desc)
    cell_nodes = CartesianCellNodes(desc.partition)
    cell_type = Fill(Int8(1),length(cell_nodes))
    new{D,T,F}(node_coords,cell_nodes,cell_type)
  end
end

OrientationStyle(::Type{<:CartesianGrid}) = Oriented()

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

function get_cell_map(grid::CartesianGrid{D,T,typeof(identity)} where {D,T})
  CartesianMap(grid.node_coords.data)
end

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

function Base.getindex(a::CartesianMap{D,T},I::Vararg{Integer,D}) where {D,T}
  p = zero(Mutable(Point{D,T}))
  x0 = a.data.origin
  dx = a.data.sizes
  @inbounds for d in 1:D
    p[d] =  x0[d] + (I[d]-1)*dx[d]
  end
  origin =  Point(p)
  grad = diagonal_tensor(VectorValue(dx))
  AffineMap(grad,origin)
end

function lazy_map(::typeof(∇),a::CartesianMap)
  dx = a.data.sizes
  grad = diagonal_tensor(VectorValue(dx))
  j = ConstantField(grad)
  Fill(j,size(a))
end

function lazy_map(::typeof(∇),a::LazyArray{<:Fill{<:Reindex{<:CartesianMap}}})
  i_to_map = a.g.value.values
  j_to_i = a.f[1]
  i_to_grad = lazy_map(∇,i_to_map)
  lazy_map(Reindex(i_to_grad),j_to_i)
end

