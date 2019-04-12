
struct CellValueFromArray{T,N,V<:AbstractArray{T,N}} <: IndexCellValue{T,N}
  v::V
end

@propagate_inbounds function getindex(self::CellValueFromArray,cell::Int)
  @inbounds self.v[cell]
end

@propagate_inbounds function getindex(self::CellValueFromArray{T,N},I::Vararg{Int,N}) where {T,N}
  @inbounds self.v[I...]
end

size(self::CellValueFromArray) = size(self.v)

IndexStyle(::Type{CellValueFromArray{T,N,V}}) where {T,N,V} = IndexStyle(V)

struct CellArrayFromArrayOfArrays{T,N,D,A,C<:AbstractArray{A,D}} <: IndexCellArray{T,N,A,D}
  c::C
end

function CellArrayFromArrayOfArrays(c::AbstractArray{A,D}) where {A<:AbstractArray,D}
  T = eltype(A)
  N = ndims(A)
  C = typeof(c)
  CellArrayFromArrayOfArrays{T,N,D,A,C}(c)
end

@propagate_inbounds function getindex(self::CellArrayFromArrayOfArrays,cell::Int)
  @inbounds self.c[cell]
end

@propagate_inbounds function getindex(self::CellArrayFromArrayOfArrays{T,N,D},I::Vararg{Int,D}) where {T,N,D}
  @inbounds self.c[I...]
end

size(self::CellArrayFromArrayOfArrays) = size(self.c)

IndexStyle(::Type{CellArrayFromArrayOfArrays{T,N,D,A,C}}) where {T,N,D,A,C} = IndexStyle(C)

mutable struct CachedSubVector{T,V<:AbstractArray{T,1}} <: AbstractArray{T,1}
  vector::V
  pini::Int
  pend::Int
end

size(self::CachedSubVector) = (1+self.pend-self.pini,)

@propagate_inbounds function getindex(self::CachedSubVector, i::Int)
  @inbounds self.vector[self.pini+i-1]
end

IndexStyle(::Type{CachedSubVector{T,V}}) where {T,V} = IndexLinear()

function locate!(self::CachedSubVector,pini,pend)
  self.pini = pini
  self.pend = pend
end

struct CellVectorFromDataAndPtrs{T,V,P} <: IndexCellArray{T,1,CachedSubVector{T,V},1}
  data::V
  ptrs::P
  cv::CachedSubVector{T,V}
end

function CellVectorFromDataAndPtrs(data::AbstractArray{T,1},ptrs::AbstractArray{Int,1}) where T
  V = typeof(data)
  P = typeof(ptrs)
  cv = CachedSubVector(data,0,0)
  CellVectorFromDataAndPtrs{T,V,P}(data,ptrs,cv)
end

@propagate_inbounds function getindex(self::CellVectorFromDataAndPtrs,cell::Int)
  pini = self.ptrs[cell]
  pend = self.ptrs[cell+1]-1
  locate!(self.cv,pini,pend)
  self.cv
end

size(self::CellVectorFromDataAndPtrs) = (length(self.ptrs)-1,)

IndexStyle(::Type{CellVectorFromDataAndPtrs{T,V,P}}) where {T,V,P} = IndexLinear()

function cellsize(self::CellVectorFromDataAndPtrs)
  l = 0
  for i in 1:(length(self.ptrs)-1)
    l = max(l,self.ptrs[i+1]-self.ptrs[i])
  end
  (l,)
end

struct CellVectorFromDataAndStride{T,V} <: IndexCellArray{T,1,CachedSubVector{T,V},1}
  data::V
  stride::Int
  cv::CachedSubVector{T,V}
end

function CellVectorFromDataAndStride(data::AbstractArray{T,1},stride::Int) where T
  V = typeof(data)
  cv = CachedSubVector(data,0,0)
  CellVectorFromDataAndStride{T,V}(data,stride,cv)
end

@propagate_inbounds function getindex(self::CellVectorFromDataAndStride,cell::Int)
  pini = self.stride*(cell - 1) + 1
  pend = self.stride*cell
  locate!(self.cv,pini,pend)
  self.cv
end

size(self::CellVectorFromDataAndStride) = (ceil(Int,length(self.data)/self.stride),)

IndexStyle(::Type{CellVectorFromDataAndStride{T,V}}) where {T,V} = IndexLinear()

function cellsize(self::CellVectorFromDataAndStride)
  (self.stride,)
end

struct CellVectorFromLocalToGlobal{T,L<:IndexCellArray{Int,1},V<:IndexCellValue{T}} <: IndexCellArray{T,1,CachedArray{T,1,Array{T,1}},1}
  lid_to_gid::L
  gid_to_val::V
  cv::CachedArray{T,1,Array{T,1}}
end

function CellVectorFromLocalToGlobal(lid_to_gid::IndexCellArray{Int,1}, gid_to_val::IndexCellValue{T}) where T
  L = typeof(lid_to_gid)
  V = typeof(gid_to_val)
  a = Vector{T}(undef,(celllength(lid_to_gid),))
  cv = CachedArray(a,size(a))
  CellVectorFromLocalToGlobal{T,L,V}(lid_to_gid,gid_to_val,cv)
end

function CellVectorFromLocalToGlobal(lid_to_gid::IndexCellArray{Int,1},gid_to_val::AbstractArray{T,N}) where {T,N}
  _gid_to_val = CellValueFromArray(gid_to_val)
  CellVectorFromLocalToGlobal(lid_to_gid,_gid_to_val)
end

@propagate_inbounds function getindex(self::CellVectorFromLocalToGlobal,cell::Int)
  lid_to_gid = self.lid_to_gid[cell]
  setsize!(self.cv,(length(lid_to_gid),))
  for (lid,gid) in enumerate(lid_to_gid)
    self.cv[lid] = self.gid_to_val[gid]
  end
  self.cv
end

size(self::CellVectorFromLocalToGlobal) = (length(self.lid_to_gid),)

IndexStyle(::Type{CellVectorFromLocalToGlobal{T,L,V}}) where {T,L,V} = IndexLinear()

cellsize(self::CellVectorFromLocalToGlobal) = cellsize(self.lid_to_gid)


