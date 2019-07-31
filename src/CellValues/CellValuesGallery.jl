module CellValuesGallery

using Gridap
using Gridap.Helpers
using Gridap.CachedSubVectors
using Gridap.CachedArrays

using Base: @propagate_inbounds, @pure

export CellValueFromArray
export CellArrayFromArrayOfArrays
export CellVectorFromDataAndPtrs
export CellVectorFromDataAndStride
export CellVectorFromLocalToGlobal
export CellVectorFromLocalToGlobalPosAndNeg
export CellVectorByComposition
import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle
import Gridap: compress

struct CellValueFromArray{T,N,V<:AbstractArray{T,N}} <: IndexCellValue{T,N}
  v::V
end

@propagate_inbounds function getindex(
  self::CellValueFromArray{T,N},I::Vararg{Int,N}) where {T,N}
  @inbounds self.v[I...]
end

size(self::CellValueFromArray) = size(self.v)

@pure IndexStyle(::Type{CellValueFromArray{T,N,V}}) where {T,N,V} = IndexStyle(V)

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

function compress(cv::CellVectorFromDataAndPtrs)
  (cv.data, cv.ptrs)
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

struct CellVectorFromLocalToGlobal{
  T,D,L<:IndexCellArray{Int,1},V<:IndexCellValue{T}} <: IndexCellArray{T,1,CachedArray{T,1,Array{T,1}},D}
  lid_to_gid::L
  gid_to_val::V
  cv::CachedArray{T,1,Array{T,1}}
end

function CellVectorFromLocalToGlobal(
  lid_to_gid::IndexCellArray{Int,1,A,D},
  gid_to_val::IndexCellValue{T}) where {T,A,D}
  L = typeof(lid_to_gid)
  V = typeof(gid_to_val)
  cv = CachedArray(T,1)
  CellVectorFromLocalToGlobal{T,D,L,V}(lid_to_gid,gid_to_val,cv)
end

function CellVectorFromLocalToGlobal(
  lid_to_gid::IndexCellArray{Int,1,A,D},
  gid_to_val::AbstractArray{T}) where {T,A,D}
  _gid_to_val = CellValueFromArray(gid_to_val)
  CellVectorFromLocalToGlobal(lid_to_gid,_gid_to_val)
end

@propagate_inbounds function getindex(self::CellVectorFromLocalToGlobal,cell::Vararg{Int,N} where N)
  lid_to_gid = self.lid_to_gid[cell...]
  setsize!(self.cv,(length(lid_to_gid),))
  for (lid,gid) in enumerate(lid_to_gid)
    self.cv[lid] = self.gid_to_val[gid]
  end
  self.cv
end

size(self::CellVectorFromLocalToGlobal) = size(self.lid_to_gid)

IndexStyle(::Type{CellVectorFromLocalToGlobal{T,D,L,V}}) where {T,D,L,V} = IndexStyle(L)

struct CellVectorFromLocalToGlobalPosAndNeg{
  T,D,L<:IndexCellArray{Int,1},V<:IndexCellValue{T},W<:IndexCellValue{T}} <: IndexCellArray{T,1,CachedArray{T,1,Array{T,1}},D}
  lid_to_gid::L
  gid_to_val_pos::V
  gid_to_val_neg::W
  cv::CachedArray{T,1,Array{T,1}}
end

function CellVectorFromLocalToGlobalPosAndNeg(
  lid_to_gid::IndexCellArray{Int,1,A,D},
  gid_to_val_pos::IndexCellValue{T},
  gid_to_val_neg::IndexCellValue{T}) where {T,A,D}
  L = typeof(lid_to_gid)
  V = typeof(gid_to_val_pos)
  W = typeof(gid_to_val_neg)
  cv = CachedArray(T,1)
  CellVectorFromLocalToGlobalPosAndNeg{T,D,L,V,W}(lid_to_gid,gid_to_val_pos,gid_to_val_neg,cv)
end

function CellVectorFromLocalToGlobalPosAndNeg(
  lid_to_gid::IndexCellArray{Int,1,A,D},
  gid_to_val_pos::AbstractArray{T},
  gid_to_val_neg::AbstractArray{T}) where {T,A,D}
  _gid_to_val_pos = CellValueFromArray(gid_to_val_pos)
  _gid_to_val_neg = CellValueFromArray(gid_to_val_neg)
  CellVectorFromLocalToGlobalPosAndNeg(lid_to_gid,_gid_to_val_pos,_gid_to_val_neg)
end

IndexStyle(::Type{CellVectorFromLocalToGlobalPosAndNeg{T,D,L,V,W}}) where {T,D,L,V,W} = IndexStyle(L)

@propagate_inbounds function getindex(
  self::CellVectorFromLocalToGlobalPosAndNeg,cell::Vararg{Int,N} where N)
  lid_to_gid = self.lid_to_gid[cell...]
  setsize!(self.cv,(length(lid_to_gid),))
  for (lid,gid) in enumerate(lid_to_gid)
    if gid > 0
      self.cv[lid] = self.gid_to_val_pos[gid]
    elseif gid < 0
      self.cv[lid] = self.gid_to_val_neg[-gid]
    else
      @unreachable
    end
  end
  self.cv
end

size(self::CellVectorFromLocalToGlobalPosAndNeg) = size(self.lid_to_gid)

struct CellVectorByComposition{
  T,L<:IndexCellVector{Int},V<:IndexCellVector{T}} <: IndexCellVector{T,CachedVector{T,Vector{T}},1}

  cell_to_x::L
  x_to_vals::V
  cv::CachedVector{T,Vector{T}}
end

function CellVectorByComposition(
  cell_to_x::IndexCellVector{Int}, x_to_vals::IndexCellVector{T}) where T
  L = typeof(cell_to_x)
  V = typeof(x_to_vals)
  cv = CachedArray(T,1)
  CellVectorByComposition{T,L,V}(cell_to_x, x_to_vals, cv)
end

@propagate_inbounds function getindex(self::CellVectorByComposition,cell::Int)
  cell_to_x = self.cell_to_x[cell]
  l = 0
  for x in cell_to_x
    l += length(self.x_to_vals[x])
  end
  setsize!(self.cv,(l,))
  l = 1
  for x in cell_to_x
    for val in self.x_to_vals[x]
      for iv in val
        self.cv[l] = iv
        l += 1
      end
    end
  end
  self.cv
end

size(self::CellVectorByComposition) = (length(self.cell_to_x),)

IndexStyle(::Type{CellVectorByComposition{T,L,V}}) where {T,L,V} = IndexLinear()

end # module
