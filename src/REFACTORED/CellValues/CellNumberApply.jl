module CellNumberApply

using Gridap
using Gridap.Helpers

export apply
import Base: iterate
import Base: length
import Base: size
import Base: getindex
import Base: IndexStyle

function apply(k::NumberKernel,v::Vararg{<:CellValue})
  CellNumberFromKernel(k,v...)
end

function apply(k::NumberKernel,v::Vararg{<:IndexCellValue})
  IndexCellNumberFromKernel(k,v...)
end

function apply(f::Function,v::Vararg;broadcast=false)
  _apply(f,v,Val(broadcast))
end

function _apply(f,v,::Val{false})
  t = _eltypes(v)
  T = Base._return_type(f,t)
  @assert T <: NumberLike
  k = NumberKernelFromFunction(f)
  apply(k,v...)
end

function _eltypes(v)
  tuple([ eltype(vi) for vi in v ]...)
end

struct CellNumberFromKernel{T,K,V} <: IterCellNumber{T}
  kernel::K
  cellvalues::V
end

function CellNumberFromKernel(k::NumberKernel,v::Vararg{<:CellValue})
  _checks(v)
  T = _compute_type(k,v)
  K = typeof(k)
  _v = Container(v...)
  V = typeof(_v)
  CellNumberFromKernel{T,K,V}(k,_v)
end

function _checks(v)
  @assert length(v) > 0
  v1, = v
  l = length(v1)
  @assert all( ( length(vi) == l for vi in v ) )
end

function _compute_type(k,v)
  t = _eltypes(v)
  T = compute_type(k,t...)
  @assert T <: NumberLike
  T
end

function length(self::CellNumberFromKernel)
  vi, = _cv(self.cellvalues)
  length(vi)
end

@inline function iterate(self::CellNumberFromKernel)
  zipped = zip(_cv(self.cellvalues)...)
  znext = iterate(zipped)
  _iterate(self,znext,zipped)
end

@inline function iterate(self::CellNumberFromKernel,state)
  zipped, zstate = state
  znext = iterate(zipped,zstate)
  _iterate(self,znext,zipped)
end

@inline function _iterate(self::CellNumberFromKernel,znext,zipped)
  if znext === nothing; return nothing end
  a, zstate = znext
  r = compute_value(self.kernel,a...)
  state = (zipped,zstate)
  (r, state)
end

struct IndexCellNumberFromKernel{T,K,V} <: IndexCellNumber{T,1}
  kernel::K
  cellvalues::V
end

function IndexCellNumberFromKernel(k::NumberKernel,v::Vararg{<:IndexCellValue})
  _checks(v)
  T = _compute_type(k,v)
  K = typeof(k)
  _v = Container(v...)
  V = typeof(_v)
  IndexCellNumberFromKernel{T,K,V}(k,_v)
end

function length(self::IndexCellNumberFromKernel)
  vi, = _cv(self.cellvalues)
  length(vi)
end

size(self::IndexCellNumberFromKernel) = (length(self),)

function getindex(self::IndexCellNumberFromKernel,i::Integer)
  vals = _getvalues(i,_cv(self.cellvalues)...)
  compute_value(self.kernel,vals...)
end

# TODO use a generated function?
_getvalues(i,v...) = @notimplemented
_getvalues(i,v1) = (v1[i],)
_getvalues(i,v1,v2) = (v1[i],v2[i])
_getvalues(i,v1,v2,v3) = (v1[i],v2[i],v3[i])
_getvalues(i,v1,v2,v3,v4) = (v1[i],v2[i],v3[i],v4[i])
_getvalues(i,v1,v2,v3,v4,v5) = (v1[i],v2[i],v3[i],v4[i],v5[i])
_getvalues(i,v1,v2,v3,v4,v5,v6) = (v1[i],v2[i],v3[i],v4[i],v5[i],v6[i])

# TODO

mutable struct Container1{V1}
  v1::V1
end

mutable struct Container2{V1,V2}
  v1::V1
  v2::V2
end

mutable struct Container3{V1,V2,V3}
  v1::V1
  v2::V2
  v3::V3
end

Container(v1) = Container1(v1)

Container(v1,v2) = Container2(v1,v2)

Container(v1,v2,v3) = Container3(v1,v2,v3)

_cv(s::Container1) = (s.v1,)

_cv(s::Container2) = (s.v1,s.v2)

_cv(s::Container3) = (s.v1,s.v2,s.v3)

function _cv!(s::Container1,v1)
  s.v1 = v1
end

function _cv!(s::Container2,v1,v2)
  s.v1 = v1
  s.v2 = v2
end

function _cv!(s::Container3,v1,v2,v3)
  s.v1 = v1
  s.v2 = v2
  s.v3 = v3
end

_container_type(V1) = Container1{V1}

_container_type(V1,V2) = Container2{V1,V2}

_container_type(V1,V2,V3) = Container3{V1,V2,V3}


end # module CellNumberApply
