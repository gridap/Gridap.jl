module CellArrayApply

using Gridap
using Gridap.Helpers
using Gridap.CachedArrays
using Gridap.CellNumberApply: _checks
using Gridap.CellNumberApply: _getvalues
using Gridap.Kernels: _size_for_broadcast

import Gridap: apply
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle
import Gridap.CellNumberApply: _apply

function apply(k::ArrayKernel,v::Vararg{<:CellValue})
  CellArrayFromKernel(k,v...)
end

function apply(k::ArrayKernel,v::Vararg{<:IndexCellValue})
  IndexCellArrayFromKernel(k,v...)
end

function _apply(f,v,::Val{true})
  k = ArrayKernelFromBroadcastedFunction(f)
  apply(k,v...)
end

struct CellArrayFromKernel{T,N,K,V} <: IterCellArray{T,N,CachedArray{T,N,Array{T,N}}}
  kernel::K
  cellvalues::V
end

function CellArrayFromKernel(k::ArrayKernel,v::Vararg{<:CellValue})
  _checks(v)
  T = _compute_T(k,v)
  N = _compute_N(k,v)
  K = typeof(k)
  V = typeof(v)
  CellArrayFromKernel{T,N,K,V}(k,v)
end

function _compute_T(k,v)
  t = _compute_eltype(v...)
  T = compute_type(k,t...)
  @assert T <: NumberLike
  T
end

function _compute_N(k,v)
  d = _compute_ndims(v...)
  compute_ndim(k,d...)
end

_nd(v::CellNumber) = 0

_nd(v::CellArray{T,N}) where {T,N} = N

# TODO use a generated function here
_compute_ndims(v...) = @notimplemented
_compute_ndims(v1) = (_nd(v1),)
_compute_ndims(v1,v2) = (_nd(v1),_nd(v2))
_compute_ndims(v1,v2,v3) = (_nd(v1),_nd(v2),_nd(v3))
_compute_ndims(v1,v2,v3,v4) = (_nd(v1),_nd(v2),_nd(v3),_nd(v4))
_compute_ndims(v1,v2,v3,v4,v5) = (_nd(v1),_nd(v2),_nd(v3),_nd(v4),_nd(v5))
_compute_ndims(v1,v2,v3,v4,v5,v6) = (_nd(v1),_nd(v2),_nd(v3),_nd(v4),_nd(v5),_nd(v6))

_eltype(v::CellNumber{T}) where T = T

_eltype(v::CellArray{T}) where T = T

const _et = _eltype

# TODO use a generated function here
_compute_eltype(v...) = @notimplemented
_compute_eltype(v1) = (_et(v1),)
_compute_eltype(v1,v2) = (_et(v1),_et(v2))
_compute_eltype(v1,v2,v3) = (_et(v1),_et(v2),_et(v3))
_compute_eltype(v1,v2,v3,v4) = (_et(v1),_et(v2),_et(v3),_et(v4))
_compute_eltype(v1,v2,v3,v4,v5) = (_et(v1),_et(v2),_et(v3),_et(v4),_et(v5))
_compute_eltype(v1,v2,v3,v4,v5,v6) = (_et(v1),_et(v2),_et(v3),_et(v4),_et(v5),_et(v6))

function length(self::CellArrayFromKernel)
  vi, = self.cellvalues
  length(vi)
end

@inline function iterate(self::CellArrayFromKernel{T,N}) where {T,N}
  zipped = zip(self.cellvalues...)
  znext = iterate(zipped)
  v = CachedArray(T,N)
  _iterate(self,znext,zipped,v)
end

@inline function iterate(self::CellArrayFromKernel,state)
  zipped, zstate, v = state
  znext = iterate(zipped,zstate)
  _iterate(self,znext,zipped,v)
end

@inline function _iterate(self::CellArrayFromKernel,znext,zipped,v)
  if znext === nothing; return nothing end
  a, zstate = znext
  s = _compute_sizes(a...)
  z = compute_size(self.kernel,s...)
  setsize!(v,z)
  compute_value!(v,self.kernel,a...)
  state = (zipped,zstate,v)
  (v, state)
end

const _bs = _size_for_broadcast

# TODO use a generated function?
_compute_sizes(a...) = @notimplemented
_compute_sizes(a1) = (_bs(a1),)
_compute_sizes(a1,a2) = (_bs(a1),_bs(a2))
_compute_sizes(a1,a2,a3) = (_bs(a1),_bs(a2),_bs(a3))
_compute_sizes(a1,a2,a3,a4) = (_bs(a1),_bs(a2),_bs(a3),_bs(a4))
_compute_sizes(a1,a2,a3,a4,a5) = (_bs(a1),_bs(a2),_bs(a3),_bs(a4),_bs(a5))
_compute_sizes(a1,a2,a3,a4,a5,a6) = (_bs(a1),_bs(a2),_bs(a3),_bs(a4),_bs(a5),_bs(a6))

struct IndexCellArrayFromKernel{T,N,K,V} <: IndexCellArray{T,N,CachedArray{T,N,Array{T,N}},1}
  kernel::K
  cellvalues::V
  cache::CachedArray{T,N,Array{T,N}}
end

function IndexCellArrayFromKernel(k::ArrayKernel,v::Vararg{<:CellValue})
  _checks(v)
  T = _compute_T(k,v)
  N = _compute_N(k,v)
  K = typeof(k)
  V = typeof(v)
  cache = CachedArray(T,N)
  IndexCellArrayFromKernel{T,N,K,V}(k,v,cache)
end

function length(self::IndexCellArrayFromKernel)
  vi, = self.cellvalues
  length(vi)
end

size(self::IndexCellArrayFromKernel) = (length(self),)

function getindex(self::IndexCellArrayFromKernel,i::Integer)
  a = _getvalues(i,self.cellvalues...)
  s = _compute_sizes(a...)
  z = compute_size(self.kernel,s...)
  v = self.cache
  setsize!(v,z)
  compute_value!(v,self.kernel,a...)
  v
end
  
end # module CellArrayApply
