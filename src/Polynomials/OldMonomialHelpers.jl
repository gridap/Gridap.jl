
function _prepare_perms(D)
  perms = zeros(Int,D,D)
  for j in 1:D
    for d in j:D
      perms[d,j] =  d-j+1
    end
    for d in 1:(j-1)
      perms[d,j] =  d+(D-j)+1
    end
  end
  perms
end

function _evaluate_nd!(
  v::AbstractVector{V},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    Kd = Val(orders[d])
    _evaluate_1d!(Monomial,Kd,c,x,d)
  end

  o = one(T)
  k = 1

  for ci in terms

    s = o
    for d in 1:dim
      @inbounds s *= c[d,ci[d]]
    end

    k = _set_value!(v,s,k)

  end

end

function _set_value!(v::AbstractVector{V},s::T,k) where {V,T}
  ncomp = num_indep_components(V)
  z = zero(T)
  @inbounds for j in 1:ncomp
    v[k] = ntuple(i -> ifelse(i == j, s, z),Val(ncomp))
    k += 1
  end
  k
end

function _set_value!(v::AbstractVector{<:Real},s,k)
    @inbounds v[k] = s
    k+1
end

function _gradient_nd!(
  v::AbstractVector{G},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  z::AbstractVector{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    Kd = Val(orders[d])
    _derivatives_1d!(Monomial,Kd,(c,g),x,d)
  end

  o = one(T)
  k = 1

  for ci in terms

    s = z
    for i in eachindex(s)
      @inbounds s[i] = o
    end
    for q in 1:dim
      for d in 1:dim
        if d != q
          @inbounds s[q] *= c[d,ci[d]]
        else
          @inbounds s[q] *= g[d,ci[d]]
        end
      end
    end

    k = _set_gradient!(v,s,k,V)

  end

end

function _set_gradient!(
  v::AbstractVector{G},s,k,::Type{<:Real}) where G

  @inbounds v[k] = s
  k+1
end

@generated function _set_gradient!(
  v::AbstractVector{G},s,k,::Type{V}) where {V,G}
  # Git blame me for readable non-generated version

  w = zero(V)
  m = Array{String}(undef, size(G))
  N_val_dims = length(size(V))
  s_size = size(G)[1:end-N_val_dims]

  body = "T = eltype(s); z = zero(T);"
  for ci in CartesianIndices(s_size)
    id = join(Tuple(ci))
    body *= "@inbounds s$id = s[$ci];"
  end

  for j in CartesianIndices(w)
    for i in CartesianIndices(m)
      m[i] = "z"
    end
    for ci in CartesianIndices(s_size)
      id = join(Tuple(ci))
      m[ci,j] = "s$id"
    end
    body *= "@inbounds v[k] = ($(join(tuple(m...), ", ")));"
    body *= "k = k + 1;"
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k))
end

# Specialization for SymTensorValue and SymTracelessTensorValue,
# necessary as long as outer(Point, V<:AbstractSymTensorValue)::G does not
# return a tensor type that implements the appropriate symmetries of the
# gradient (and hessian)
@generated function _set_gradient!(
  v::AbstractVector{G},s,k,::Type{V}) where {V<:AbstractSymTensorValue{D},G} where D
  # Git blame me for readable non-generated version

  T = eltype(s)
  m = Array{String}(undef, size(G))
  s_length = size(G)[1]

  is_traceless = V <: SymTracelessTensorValue
  skip_last_diagval = is_traceless ? 1 : 0    # Skid V_DD if traceless

  body = "z = $(zero(T));"
  for i in 1:s_length
    body *= "@inbounds s$i = s[$i];"
  end

  for c in 1:(D-skip_last_diagval) # Go over cols
    for r in c:D                   # Go over lower triangle, current col
      for i in eachindex(m)
        m[i] = "z"
      end
      for i in 1:s_length # indices of the Vector s
        m[i,r,c] = "s$i"
        if (r!=c)
          m[i,c,r] = "s$i"
        elseif is_traceless # V_rr contributes negatively to V_DD (tracelessness)
          m[i,D,D] = "-s$i"
        end
      end
      body *= "@inbounds v[k] = ($(join(tuple(m...), ", ")));"
      body *= "k = k + 1;"
    end
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k))
end

function _hessian_nd!(
  v::AbstractVector{G},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    Kd = Val(orders[d])
    _derivatives_1d!(Monomial,Kd,(c,g,h),x,d)
  end

  z = zero(Mutable(TensorValue{D,D,T}))
  o = one(T)
  k = 1

  for ci in terms

    s = z
    for i in eachindex(s)
      @inbounds s[i] = o
    end
    for r in 1:dim
      for q in 1:dim
        for d in 1:dim
          if d != q && d != r
            @inbounds s[r,q] *= c[d,ci[d]]
          elseif d == q && d ==r
            @inbounds s[r,q] *= h[d,ci[d]]
          else
            @inbounds s[r,q] *= g[d,ci[d]]
          end
        end
      end
    end

    k = _set_gradient!(v,s,k,V)

  end

end

#function evaluate!(
#  #cache, f::FieldGradientArray{N,<:AbstractVector{Monomial}}, x::VectorValue) where {N}
#  cache, f::FieldGradientArray{N,<:AbstractVector{<:Type{Polynomial}}}, x::VectorValue) where {N}
#  r, cf, xs = cache
#  xs[1] = x
#  v = evaluate!(cf,f,xs)
#  ndof = size(v,2)
#  setsize!(r,(ndof,))
#  a = r.array
#  copyto!(a,v)
#  a
#end

# Field implementation
#function _evaluate!(
#  cache,
#  fg::FieldGradientArray{1,MonomialBasis{D,V}},
#  x::AbstractVector{<:Point},
#  TisbitsType::Val{true}) where {D,V}
#
#  f = fg.fa
#  r, v, c, g = cache
#  z = zero(Mutable(VectorValue{D,eltype(T)}))
#  np = length(x)
#  ndof = length(f)
#  n = 1 + _maximum(f.orders)
#  setsize!(r,(np,ndof))
#  setsize!(v,(ndof,))
#  setsize!(c,(D,n))
#  setsize!(g,(D,n))
#  for i in 1:np
#    @inbounds xi = x[i]
#    _gradient_nd!(v,xi,f.orders,f.terms,c,g,z,T)
#    for j in 1:ndof
#      @inbounds r[i,j] = v[j]
#    end
#  end
#  r.array
#end
#
#function _evaluate!(
#  cache,
#  fg::FieldGradientArray{1,MonomialBasis{D,V}},
#  x::AbstractVector{<:Point},
#  TisbitsType::Val{false}) where {D,V}
#
#  f = fg.fa
#  r, v, c, g, z = cache
#  np = length(x)
#  ndof = length(f)
#  n = 1 + _maximum(f.orders)
#  setsize!(r,(np,ndof))
#  setsize!(v,(ndof,))
#  setsize!(c,(D,n))
#  setsize!(g,(D,n))
#  for i in 1:np
#    @inbounds xi = x[i]
#    _gradient_nd!(v,xi,f.orders,f.terms,c,g,z,T)
#    for j in 1:ndof
#      @inbounds r[i,j] = v[j]
#    end
#  end
#  r.array
#end
#
#function evaluate!(
#  cache,
#  fg::FieldGradientArray{1,MonomialBasis{D,V}},
#  x::AbstractVector{<:Point}) where {D,V}
#
#  r, v, c, g = cache
#  TisbitsType = Val(isbitstype(eltype(c)))
#  _evaluate!(cache,fg,x,TisbitsType)
#end
#
#function return_cache(
#  fg::FieldGradientArray{2,MonomialBasis{D,V}},
#  x::AbstractVector{<:Point}) where {D,V}
#
#  f = fg.fa
#  @check D == length(eltype(x)) "Incorrect number of point components"
#  np = length(x)
#  ndof = length(f)
#  xi = testitem(x)
#  T = gradient_type(gradient_type(V,xi),xi)
#  n = 1 + _maximum(f.orders)
#  r = CachedArray(zeros(T,(np,ndof)))
#  v = CachedArray(zeros(T,(ndof,)))
#  c = CachedArray(zeros(eltype(T),(D,n)))
#  g = CachedArray(zeros(eltype(T),(D,n)))
#  h = CachedArray(zeros(eltype(T),(D,n)))
#  (r, v, c, g, h)
#end
#
#function evaluate!(
#  cache,
#  fg::FieldGradientArray{2,MonomialBasis{D,V}},
#  x::AbstractVector{<:Point}) where {D,V}
#
#  f = fg.fa
#  r, v, c, g, h = cache
#  np = length(x)
#  ndof = length(f)
#  n = 1 + _maximum(f.orders)
#  setsize!(r,(np,ndof))
#  setsize!(v,(ndof,))
#  setsize!(c,(D,n))
#  setsize!(g,(D,n))
#  setsize!(h,(D,n))
#  for i in 1:np
#    @inbounds xi = x[i]
#    _hessian_nd!(v,xi,f.orders,f.terms,c,g,h,T)
#    for j in 1:ndof
#      @inbounds r[i,j] = v[j]
#    end
#  end
#  r.array
#end

