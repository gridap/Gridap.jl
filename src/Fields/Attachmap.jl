
"""
    attachmap(f,phi)
"""
function attachmap(f,phi)
  k = MapGrad()
  apply_kernel_to_field(k,f,phi)
end

"""
    attachmap(f::AbstractArray,phi::AbstractArray)
"""
function attachmap(f::AbstractArray,phi::AbstractArray)
  k = AddMap()
  apply(k,f,phi)
end

struct MapGrad <: Kernel end

@inline apply_kernel!(cache,k::MapGrad,fx,phix) = fx

function apply_kernel_gradient(k::MapGrad,f,phi)
  g = field_gradient(f)
  jac = field_gradient(phi)
  k = PhysGrad()
  apply_kernel_to_field(k,g,jac)
end

struct PhysGrad <: Kernel end

function kernel_cache(k::PhysGrad,a,b)
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

@inline function apply_kernel!(r,k::PhysGrad,a,b)
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

struct AddMap <: Kernel end

@inline function apply_kernel!(cache,k::AddMap,fi,phii)
  attachmap(fi,phii)
end

function apply_gradient(k::AddMap,f,phi)
  g = gradient(f)
  jac = gradient(phi)
  k = PhysGrad()
  apply_to_field_array(k,g,jac)
end

function kernel_evaluate(k::AddMap,x,f,phi)
  evaluate(f,x)
end

