module MapApply

using Gridap
using Gridap.Helpers
using Gridap.CachedArrays
using Gridap.Kernels: _compute_N, _compute_T

import Gridap: apply
import Gridap: evaluate!
import Gridap: return_size

function apply(k::ArrayKernel,m::Map,v::Vararg)
  MapFromKernel(k,m,v...)
end

mutable struct MapFromKernel{S,M,T,N,K,V,C} <: Map{S,M,T,N}
  kernel::K
  inputs::V
  caches::C
end

function MapFromKernel(k::ArrayKernel,v::Vararg)
  S = _compute_S(v)
  M = _compute_M(v)
  T = _compute_T(k,v)
  N = _compute_N(k,v)
  K = typeof(k)
  V = typeof(v)
  c = _create_caches(v)
  C = typeof(c)
  MapFromKernel{S,M,T,N,K,V,C}(k,v,c)
end

function evaluate!(
  m::MapFromKernel{S,M,T,N},
  points::AbstractArray{<:S,M},
  v::AbstractArray{T,N}) where {S,M,T,N}
  _evaluate_inputs!(points,m.caches...,m.inputs...)
  compute_value!(v,m.kernel,m.caches...)
end

function return_size(m::MapFromKernel{S,M},s::NTuple{M,Int}) where {S,M}
  z = _return_sizes(s,m.inputs...)
  compute_size(m.kernel,z...)
end

function setinputs!(m::MapFromKernel,a::Vararg)
  m.inputs = a
end

function _create_caches(v)
  tuple([_cache(vi) for vi in v]...)
end

_cache(v::Map{S,M,T,N}) where {S,M,T,N} = CachedArray(T,N)

_cache(v::AbstractArray{T,N}) where {T,N} = CachedArray(T,N)

function _compute_S(v)
  @assert length(v) > 0
  v1, = v
  S = _stype(v1)
  @assert S != nothing
  @assert all([ (_stype(vi)==S || _stype(vi)==nothing ) for vi in v ])
  S
end

function _compute_M(v)
  @assert length(v) > 0
  v1, = v
  M = _m(v1)
  @assert M != nothing
  @assert all([ (_m(vi)==M || _m(vi)==nothing ) for vi in v ])
  M
end

_stype(v::Map{S}) where S = S

_stype(v::AbstractArray) = nothing

_m(v::Map{S,M}) where {S,M} = M

_m(v::AbstractArray) = nothing

_rz(s,i) = @unreachable

_rz(s,i::Map) = return_size(i,s)

_rz(s,i::AbstractArray) = size(i)

# TODO also with a generated function
_return_sizes(s,i...) = @notimplemented
_return_sizes(s,i1) = (_rz(s,i1),)
_return_sizes(s,i1,i2) = (_rz(s,i1),_rz(s,i2))
_return_sizes(s,i1,i2,i3) = (_rz(s,i1),_rz(s,i2),_rz(s,i3))
_return_sizes(s,i1,i2,i3,i4) = (_rz(s,i1),_rz(s,i2),_rz(s,i3),_rz(s,i4))
_return_sizes(s,i1,i2,i3,i4,i5) = (_rz(s,i1),_rz(s,i2),_rz(s,i3),_rz(s,i4),_rz(s,i5))
_return_sizes(s,i1,i2,i3,i4,i5,i6) = (_rz(s,i1),_rz(s,i2),_rz(s,i3),_rz(s,i4),_rz(s,i5),_rz(s,i6))

_ev!(p,c,i) = @unreachable

function _ev!(p,c,i::Map)
  s = return_size(i,size(p))
  setsize!(c,s)
  evaluate!(i,p,c)
end

function _ev!(p,c,i::AbstractArray)
  s = size(i)
  setsize!(c,s)
  for k in eachindex(i)
    @inbounds c[k] = i[k]
  end
end

# TODO also with a generated function
_evaluate_inputs!(p,ci...) = @notimplemented

function _evaluate_inputs!(p,c1,i1)
  _ev!(p,c1,i1)
end

function _evaluate_inputs!(p,c1,c2,i1,i2)
  _ev!(p,c1,i1)
  _ev!(p,c2,i2)
end

function _evaluate_inputs!(p,c1,c2,c3,i1,i2,i3)
  _ev!(p,c1,i1)
  _ev!(p,c2,i2)
  _ev!(p,c3,i3)
end

function _evaluate_inputs!(p,c1,c2,c3,c4,i1,i2,i3,i4)
  _ev!(p,c1,i1)
  _ev!(p,c2,i2)
  _ev!(p,c3,i3)
  _ev!(p,c4,i4)
end

function _evaluate_inputs!(p,c1,c2,c3,c4,c5,i1,i2,i3,i4,i5)
  _ev!(p,c1,i1)
  _ev!(p,c2,i2)
  _ev!(p,c3,i3)
  _ev!(p,c4,i4)
  _ev!(p,c5,i5)
end

function _evaluate_inputs!(p,c1,c2,c3,c4,c5,c6,i1,i2,i3,i4,i5,i6)
  _ev!(p,c1,i1)
  _ev!(p,c2,i2)
  _ev!(p,c3,i3)
  _ev!(p,c4,i4)
  _ev!(p,c5,i5)
  _ev!(p,c6,i6)
end

end # module MapApply
