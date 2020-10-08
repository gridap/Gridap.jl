
"""
    lincomb(a::Field,b::AbstractVector)

Returns a field obtained by the "linear combination" of
the value of the field basis `a` and the coefficient vector `b`.  The value of the resulting field
evaluated at a vector of points `x` is defined as

    ax = evaluate(a,x)
    ax*b

On the other hand, the gradient of the resulting field is defined as

    ∇ax = evaluate(gradient(a),x)
    ∇ax*b

"""
function lincomb(a::Field,b::AbstractVector)
  LinComField(a,b)
end

struct LinCom <: Mapping end

function return_cache(k::LinCom,a,b)
  Ta = eltype(a)
  Tb = eltype(b)
  T = return_type(outer,Ta,Tb)
  np = length(b)
  r = zeros(T,np)
  CachedArray(r)
end

function testitem!(r,k::LinCom,a,b)
  if _lincomb_valid_checks(a,b)
    evaluate!(r,k,a,b)
  else
    r.array
  end
end

function _lincomb_checks(a,b)
  s = "lincom: Number of fields in basis needs to be equal to number of coefs."
  @assert _lincomb_valid_checks(a,b) s
end

function _lincomb_valid_checks(a,b)
  nb = length(b)
  np, na = size(a)
  nb == na
end

@inline function evaluate!(r,k::LinCom,a,b)
  _lincomb_checks(a,b)
  np, nf = size(a)
  setsize!(r,(np,))
  z = zero(eltype(r))
  _r = r.array
  for i in 1:np
    @inbounds _r[i] = z
    for j in 1:nf
      @inbounds _r[i] += outer(a[i,j],b[j])
    end
  end
  _r
end

mutable struct LinComField{A,B} <: Field
  basis::A
  coefs::B
  @inline function LinComField(basis,coefs)
    A = typeof(basis)
    B = typeof(coefs)
    new{A,B}(basis,coefs)
  end
end

function return_cache(f::LinComField,x)
  ca = return_cache(f.basis,x)
  a = evaluate!(ca,f.basis,x)
  b = f.coefs
  k = LinCom()
  ck = return_cache(k,a,b)
  (ca,ck)
end

@inline function evaluate!(cache,f::LinComField,x)
  ca, ck = cache
  a = evaluate!(ca,f.basis,x)
  b = f.coefs
  k = LinCom()
  evaluate!(ck,k,a,b)
end

function field_gradient(f::LinComField)
  g = field_gradient(f.basis)
  LinComField(g,f.coefs)
end

struct LinComValued <: Mapping end

@inline function return_cache(k::LinComValued,a,b)
  LinComField(a,b)
end

@inline function evaluate!(f,k::LinComValued,a,b)
  f.basis = a
  f.coefs = b
  f
end

function lazy_map_gradient(k::LinComValued,a,b)
  g = field_array_gradient(a)
  lincomb(g,b)
end

function kernel_evaluate(k::LinComValued,x,a,b)
  ax = evaluate_field_array(a,x)
  lazy_map_lincomb(ax,b)
end

"""
    lazy_map_lincomb(ax,b)
"""
function lazy_map_lincomb(ax,b)
  k = LinCom()
  lazy_map(k,ax,b)
end

"""
    lincomb(a::AbstractArray{<:Field},b::AbstractArray)

Returns an array of fields numerically equivalent to

    map(lincomb,a,b)
"""
function lincomb(a::AbstractArray,b::AbstractArray)
  k = LinComValued()
  lazy_map(k,a,b)
end
