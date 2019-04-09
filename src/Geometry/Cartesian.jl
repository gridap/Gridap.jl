
struct CartesianGrid{D} <: Grid{D,D}
  dim_to_limits::NTuple{D,NTuple{2,Float64}}
  dim_to_ncells::NTuple{D,Int}
  extrusion::NTuple{D,Int}
end

function CartesianGrid(;domain::NTuple{D2,Float64},partition::NTuple{D,Int}) where {D2,D}
  @assert D2 == 2*D
  dim_to_limits = tuple([(domain[2*i-1],domain[2*i]) for i in 1:D ]...)
  extrusion = tuple(fill(HEX_AXIS,D)...)
  dim_to_ncells = partition
  CartesianGrid{D}(dim_to_limits,dim_to_ncells,extrusion)
end

function coordinates(self::CartesianGrid)
  dim_to_npoint = tuple([ i+1 for i in self.dim_to_ncells ]...)
  CartesianGridCoords(self.dim_to_limits,dim_to_npoint)
end

connectivity(self::CartesianGrid) = CartesianGridConnectivity(self.dim_to_ncells)

celltypes(self::CartesianGrid) = ConstantCellValue(self.extrusion,prod(self.dim_to_ncells))

# Ancillary types

struct CartesianGridCoords{D} <: IndexCellValue{Point{D},D}
  dim_to_limits::NTuple{D,NTuple{2,Float64}}
  dim_to_npoint::NTuple{D,Int}
end

size(self::CartesianGridCoords) = self.dim_to_npoint

IndexStyle(::Type{CartesianGridCoords{D}} where D) = IndexCartesian()

function getindex(self::CartesianGridCoords{D}, I::Vararg{Int, D}) where D
  p = zero(MPoint{D})
  @inbounds for d in 1:D
    xa = self.dim_to_limits[d][1]
    xb = self.dim_to_limits[d][2]
    p[d] =  xa + (I[d]-1)*(xb-xa)/(self.dim_to_npoint[d]-1)
  end
  Point{D}(p)
end

struct CartesianGridConnectivity{D,L} <: IndexCellArray{Int,1,SVector{L,Int},D}
  dim_to_ncell::SVector{D,Int}
end

function CartesianGridConnectivity(dim_to_ncell::NTuple{D,Int}) where D
  CartesianGridConnectivity{D,2^D}(dim_to_ncell)
end

size(self::CartesianGridConnectivity) = self.dim_to_ncell.data

IndexStyle(::Type{CartesianGridConnectivity{D,L}} where {D,L}) = IndexCartesian()

function getindex(self::CartesianGridConnectivity{D,L}, I::Vararg{Int, D}) where {D,L}
  dim_to_ngpoint = 1 .+ self.dim_to_ncell
  dim_to_nlpoint = @SVector fill(2,D)
  offset = @SVector fill(1,D)
  pointgids = LinearIndices(dim_to_ngpoint.data)
  cellpointlids = CartesianIndices(dim_to_nlpoint.data)
  cellgid = CartesianIndex(I...) - CartesianIndex(offset.data)
  cellpointgids = cellpointlids .+ cellgid
  ids = zero(MVector{L,Int})
  @inbounds for (l,pgid) in enumerate(cellpointgids)
    ids[l] = pointgids[pgid]
  end
  SVector{L,Int}(ids)
end

