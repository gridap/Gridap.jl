
# Descriptor of a cartesian grid

"""
    struct CartesianDescriptor{D,T,F<:Function}
      origin::Point{D,T}
      sizes::Point{D,T}
      partition::Point{D,Int}
      map::F
    end

Struct that stores the data defining a Cartesian grid.
"""
struct CartesianDescriptor{D,T,F<:Function} <: GridapType
  origin::Point{D,T}
  sizes::Point{D,T}
  partition::Point{D,Int}
  map::F
  @doc """
      CartesianDescriptor(origin,sizes,partition,map::Function=identity)
  """
  function CartesianDescriptor(origin,sizes,partition,map::Function=identity)
    D = length(partition)
    T = eltype(sizes)
    F = typeof(map)
    new{D,T,F}(origin,sizes,partition,map)
  end
end

"""
    CartesianDescriptor(domain,partition,map::Function=identity)
"""
function CartesianDescriptor(domain,partition,map::Function=identity)
  D = length(partition)
  limits = [(domain[2*d-1],domain[2*d]) for d in 1:D]
  sizes = [(limits[d][2]-limits[d][1])/partition[d] for d in 1:D]
  origin = [ limits[d][1] for d in 1:D]
  CartesianDescriptor(origin,sizes,partition,map)
end

"""
    CartesianDescriptor(
      pmin::Point{D},pmax::Point{D},partition,map::Function=identity) where D

"""
function CartesianDescriptor(
  pmin::Point{D},pmax::Point{D},partition,map::Function=identity) where D
  T = eltype(pmin)
  domain = zeros(T,2*D)
  for d in 1:D
    domain[2*(d-1)+1] = pmin[d]
    domain[2*(d-1)+2] = pmax[d]
  end
  CartesianDescriptor(domain,partition,map)
end

# Coordinates

struct CartesianCoordinates{D,T,F} <: AbstractArray{Point{D,T},D}
  data::CartesianDescriptor{D,T,F}
end

Base.size(a::CartesianCoordinates) = Tuple(a.data.partition) .+ 1

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
  partition::Point{D,Int}
  function CartesianCellNodes(partition)
    D = length(partition)
    new{D}(partition)
  end
end

Base.size(a::CartesianCellNodes) = Tuple(a.partition)

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
    CartesianGrid(domain,partition,map::Function=identity)
"""
function CartesianGrid(domain,partition,map::Function=identity)
  desc = CartesianDescriptor(domain,partition,map)
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

Base.size(a::CartesianMap) = Tuple(a.data.partition)

Base.IndexStyle(::Type{<:CartesianMap}) = IndexCartesian()

function Base.getindex(a::CartesianMap{D,T},I::Vararg{Int,D}) where {D,T}
  p = zero(mutable(Point{D,T}))
  x0 = a.data.origin
  dx = a.data.sizes
  @inbounds for d in 1:D
    p[d] =  x0[d] + (I[d]-1)*dx[d]
  end
  origin =  Point(p)
  jacobian = diagonal_tensor(dx)
  AffineMap(jacobian,origin)
end

function field_array_gradient(a::CartesianMap)
  dx = a.data.sizes
  jacobian = diagonal_tensor(dx)
  j = AffineMapGrad(jacobian)
  Fill(j,length(a))
end

function get_cell_map(grid::CartesianGrid{D,T,typeof(identity)} where {D,T})
  CartesianMap(grid.node_coords.data)
end

