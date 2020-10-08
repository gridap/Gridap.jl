
"""
    attachmap(f,phi)
"""
function attachmap(f,phi)
  k = MapGrad()
  evaluate_to_field(k,f,phi)
end

"""
    attachmap(f::AbstractArray,phi::AbstractArray)
"""
function attachmap(f::AbstractArray,phi::AbstractArray)
  k = AddMap()
  lazy_map(k,f,phi)
end

struct MapGrad <: Mapping end

@inline evaluate!(cache,k::MapGrad,fx,phix) = fx

function evaluate_gradient(k::MapGrad,f,phi)
  g = field_gradient(f)
  jac = field_gradient(phi)
  k = PhysGrad()
  evaluate_to_field(k,g,jac)
end

struct PhysGrad <: Mapping end

function return_cache(k::PhysGrad,a,b)
  _attachmap_checks(a,b)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(⋅,return_type(inv,Tb),Ta)
  r = zeros(T,size(a))
  CachedArray(r)
end

function _attachmap_checks(a,b)
  @assert ndims(a) == 2 "attachmap: map can only be attached to bases of fields."
  ia,ja = size(a)
  ib = length(b)
  @assert ia == ib "attachmap: basis and jacobian size mismatch."
end

@inline function evaluate!(r,k::PhysGrad,a,b)
  _attachmap_checks(a,b)
  s = size(a)
  setsize!(r,s)
  c = r.array
  np, ni = s
  for p in 1:np
    @inbounds jacinv = inv(b[p])
    for i in 1:ni
      @inbounds c[p,i] = jacinv ⋅ a[p,i]
    end
  end
  c
end

struct AddMap <: Mapping end

@inline function evaluate!(cache,k::AddMap,fi,phii)
  attachmap(fi,phii)
end

function lazy_map_gradient(k::AddMap,f,phi)
  g = gradient(f)
  jac = gradient(phi)
  k = PhysGrad()
  lazy_map_to_field_array(k,g,jac)
end

function kernel_evaluate(k::AddMap,x,f,phi)
  evaluate(f,x)
end
