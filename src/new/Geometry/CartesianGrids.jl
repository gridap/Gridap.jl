
# Descriptor of a cartesian grid

struct CartesianDescriptor{D,T,F<:Function}
  origin::Point{D,T}
  sizes::Point{D,T}
  partition::Point{D,Int}
  map::F
  function CartesianDescriptor(origin,sizes,partition,map::Function=identity)
    D = length(partition)
    T = eltype(sizes)
    F = typeof(map)
    new{D,T,F}(origin,sizes,partition,map)
  end
end

function CartesianDescriptor(domain,partition,map::Function=identity)
  D = length(partition)
  limits = [(domain[2*d-1],domain[2*d]) for d in 1:D]
  sizes = [(limits[d][2]-limits[d][1])/partition[d] for d in 1:D]
  origin = [ limits[d][1] for d in 1:D]
  CartesianDescriptor(origin,sizes,partition,map)
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
"""
struct CartesianGrid{D,T,F} <: ConformingTriangulation{D,D}
  node_coords::CartesianCoordinates{D,T,F}
  cell_nodes::CartesianCellNodes{D}
  cell_type::Fill{Int8,1,Tuple{Base.OneTo{Int}}}
  function CartesianGrid(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
    node_coords = CartesianCoordinates(desc)
    cell_nodes = CartesianCellNodes(desc.partition)
    cell_type = Fill(Int8(1),length(cell_nodes))
    new{D,T,F}(node_coords,cell_nodes,cell_type)
  end
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

"""
"""
function CartesianGrid(domain,partition,map=identity)
  desc = CartesianDescriptor(domain,partition,map)
  CartesianGrid(desc)
end


