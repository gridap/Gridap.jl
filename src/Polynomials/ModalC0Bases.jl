struct ModalC0BasisFunction <: Field end

struct ModalC0Basis{D,T,V} <: AbstractVector{ModalC0BasisFunction}
  orders::NTuple{D,Int}
  terms::Vector{CartesianIndex{D}}
  a::Vector{Point{D,V}}
  b::Vector{Point{D,V}}
  function ModalC0Basis{D}(
    ::Type{T},
    orders::NTuple{D,Int},
    terms::Vector{CartesianIndex{D}},
    a::Vector{Point{D,V}},
    b::Vector{Point{D,V}}) where {D,T,V}
    new{D,T,V}(orders,terms,a,b)
  end
end

@inline Base.size(a::ModalC0Basis{D,T,V}) where {D,T,V} = (length(a.terms)*num_components(T),)
@inline Base.getindex(a::ModalC0Basis,i::Integer) = ModalC0BasisFunction()
@inline Base.IndexStyle(::ModalC0Basis) = IndexLinear()

function ModalC0Basis{D}(
  ::Type{T},
  orders::NTuple{D,Int},
  a::Vector{Point{D,V}},
  b::Vector{Point{D,V}};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,T,V}

  terms = _define_terms_mc0(filter, sort!, orders)
  ModalC0Basis{D}(T,orders,terms,a,b)
end

function ModalC0Basis{D}(
  ::Type{T},
  orders::NTuple{D,Int},
  sa::Point{D,V},
  sb::Point{D,V};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,T,V}

  terms = _define_terms_mc0(filter, sort!, orders)
  a = fill(sa,length(terms))
  b = fill(sb,length(terms))
  ModalC0Basis{D}(T,orders,terms,a,b)
end

function ModalC0Basis{D}(
  ::Type{T},
  orders::NTuple{D,Int};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,T}

  sa = Point{D,eltype(T)}(tfill(zero(eltype(T)),Val{D}()))
  sb = Point{D,eltype(T)}(tfill(one(eltype(T)),Val{D}()))
  ModalC0Basis{D}(T,orders,sa,sb,filter=filter,sort! = sort!)
end

function ModalC0Basis{D}(
  ::Type{T},
  order::Int,
  a::Vector{Point{D,V}},
  b::Vector{Point{D,V}};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,T,V}

  orders = tfill(order,Val{D}())
  ModalC0Basis{D}(T,orders,a,b,filter=filter,sort! = sort!)
end

function ModalC0Basis{D}(
  ::Type{T},
  order::Int;
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,T}

  orders = tfill(order,Val{D}())
  ModalC0Basis{D}(T,orders,filter=filter,sort! = sort!)
end

# API

"""
    get_order(b::ModalC0Basis)
"""
function get_order(b::ModalC0Basis)
  maximum(b.orders)
end

"""
    get_orders(b::ModalC0Basis)
"""
function get_orders(b::ModalC0Basis)
  b.orders
end

return_type(::ModalC0Basis{D,T,V}) where {D,T,V} = T

# Field implementation

function return_cache(f::ModalC0Basis{D,T,V},x::AbstractVector{<:Point}) where {D,T,V}
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f.terms)*num_components(T)
  n = 1 + _maximum(f.orders)
  r = CachedArray(zeros(T,(np,ndof)))
  v = CachedArray(zeros(T,(ndof,)))
  c = CachedArray(zeros(eltype(T),(D,n)))
  (r, v, c)
end

function evaluate!(cache,f::ModalC0Basis{D,T,V},x::AbstractVector{<:Point}) where {D,T,V}
  r, v, c = cache
  np = length(x)
  ndof = length(f.terms)*num_components(T)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _evaluate_nd_mc0!(v,xi,f.a,f.b,f.orders,f.terms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{1,ModalC0Basis{D,V,W}},
  x::AbstractVector{<:Point}) where {D,V,W}

  f = fg.fa
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f.terms)*num_components(V)
  xi = testitem(x)
  T = gradient_type(V,xi)
  n = 1 + _maximum(f.orders)
  r = CachedArray(zeros(T,(np,ndof)))
  v = CachedArray(zeros(T,(ndof,)))
  c = CachedArray(zeros(eltype(T),(D,n)))
  g = CachedArray(zeros(eltype(T),(D,n)))
  (r, v, c, g)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{1,ModalC0Basis{D,T,V}},
  x::AbstractVector{<:Point}) where {D,T,V}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = length(f.terms) * num_components(T)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_mc0!(v,xi,f.a,f.b,f.orders,f.terms,c,g,T)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{2,ModalC0Basis{D,V,W}},
  x::AbstractVector{<:Point}) where {D,V,W}

  f = fg.fa
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = length(f.terms)*num_components(V)
  xi = testitem(x)
  T = gradient_type(gradient_type(V,xi),xi)
  n = 1 + _maximum(f.orders)
  r = CachedArray(zeros(T,(np,ndof)))
  v = CachedArray(zeros(T,(ndof,)))
  c = CachedArray(zeros(eltype(T),(D,n)))
  g = CachedArray(zeros(eltype(T),(D,n)))
  h = CachedArray(zeros(eltype(T),(D,n)))
  (r, v, c, g, h)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{2,ModalC0Basis{D,T,V}},
  x::AbstractVector{<:Point}) where {D,T,V}

  f = fg.fa
  r, v, c, g, h = cache
  np = length(x)
  ndof = length(f.terms) * num_components(T)
  n = 1 + _maximum(f.orders)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  setsize!(h,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
    _hessian_nd_mc0!(v,xi,f.a,f.b,f.orders,f.terms,c,g,h,T)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

# Optimizing evaluation at a single point

function return_cache(f::AbstractVector{ModalC0BasisFunction},x::Point)
  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(cache,f::AbstractVector{ModalC0BasisFunction},x::Point)
  r, cf, xs = cache
  xs[1] = x
  v = evaluate!(cf,f,xs)
  ndof = size(v,2)
  setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

function return_cache(
  f::FieldGradientArray{N,<:AbstractVector{ModalC0BasisFunction}}, x::Point) where {N}
  xs = [x]
  cf = return_cache(f,xs)
  v = evaluate!(cf,f,xs)
  r = CachedArray(zeros(eltype(v),(size(v,2),)))
  r, cf, xs
end

function evaluate!(
  cache, f::FieldGradientArray{N,<:AbstractVector{ModalC0BasisFunction}}, x::Point) where {N}
  r, cf, xs = cache
  xs[1] = x
  v = evaluate!(cf,f,xs)
  ndof = size(v,2)
  setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

# Helpers

_s_filter_mc0(e,o) = ( sum( [ i for i in e if i>1 ] ) <= o )

_sort_by_tensor_prod!(terms,orders) = terms

function _sort_by_nfaces!(terms::Vector{CartesianIndex{D}},orders) where D

  # Generate indices of n-faces and order s.t.
  # (1) dimension-increasing (2) lexicographic
  bin_rang_nfaces = tfill(0:1,Val{D}())
  bin_ids_nfaces = vec(collect(Iterators.product(bin_rang_nfaces...)))
  sum_bin_ids_nfaces = sum.(bin_ids_nfaces)
  bin_ids_nfaces = permute!(bin_ids_nfaces,sortperm(sum_bin_ids_nfaces))

  # Generate LIs of basis funs s.t. order by n-faces
  lids_b = LinearIndices(Tuple([orders[i]+1 for i=1:D]))

  eet = eltype(eltype(bin_ids_nfaces))
  f(x) = Tuple( x[i] == one(eet) ? (0:0) : (1:2) for i in 1:length(x) )
  g(x) = Tuple( x[i] == one(eet) ? (3:orders[i]+1) : (0:0) for i in 1:length(x) )
  rang_nfaces = map(f,bin_ids_nfaces)
  rang_own_dofs = map(g,bin_ids_nfaces)

  P = Int64[]
  for i = 1:length(bin_ids_nfaces)
    cis_nfaces = CartesianIndices(rang_nfaces[i])
    cis_own_dofs = CartesianIndices(rang_own_dofs[i])
    for ci in cis_nfaces
      ci = ci .+ cis_own_dofs
      P = vcat(P,reshape(lids_b[ci],length(ci)))
    end
  end

  permute!(terms,P)
end

function _compute_filter_mask(terms,filter,orders)
  g = (0 .* orders) .+ 1
  to = CartesianIndex(g)
  maxorder = _maximum(orders)
  term_to_is_fterm = lazy_map(t->filter(Int[Tuple(t-to)...],maxorder),terms)
  findall(term_to_is_fterm)
end

function _define_terms_mc0(filter,sort!,orders)
  terms = _define_terms(_q_filter,orders)
  sort!(terms,orders)
  mask = _compute_filter_mask(terms,filter,orders)
  collect(lazy_map(Reindex(terms),mask))
end

function _evaluate_1d_mc0!(v::AbstractMatrix{T},x,a,b,order,d) where T
  @assert order > 0
  n = order + 1
  z = one(T)
  @inbounds v[d,1] = z - x[d]
  @inbounds v[d,2] = x[d]
  if n > 2
    ξ = ( 2*x[d] - ( a[d] + b[d] ) ) / ( b[d] - a[d] )
    for i in 3:n
      @inbounds v[d,i] = -sqrt(2*i-3)*v[d,1]*v[d,2]*jacobi(ξ,i-3,1,1)/(i-2)
    end
  end
end

function _gradient_1d_mc0!(v::AbstractMatrix{T},x,a,b,order,d) where T
  @assert order > 0
  n = order + 1
  z = one(T)
  @inbounds v[d,1] = -z
  @inbounds v[d,2] = z
  if n > 2
    ξ = ( 2*x[d] - ( a[d] + b[d] ) ) / ( b[d] - a[d] )
    v1 = z - x[d]
    v2 = x[d]
    for i in 3:n
      j, dj = jacobi_and_derivative(ξ,i-3,1,1)
      @inbounds v[d,i] = -sqrt(2*i-3)*(v[d,1]*v2*j+v1*v[d,2]*j+v1*v2*(2/(b[d]-a[d]))*dj)/(i-2)
    end
  end
end

function _hessian_1d_mc0!(v::AbstractMatrix{T},x,a,b,order,d) where T
  @assert order > 0
  n = order + 1
  y = zero(T)
  z = one(T)
  @inbounds v[d,1] = y
  @inbounds v[d,2] = y
  if n > 2
    ξ = ( 2*x[d] - ( a[d] + b[d] ) ) / ( b[d] - a[d] )
    v1 = z - x[d]
    v2 = x[d]
    dv1 = -z
    dv2 = z
    for i in 3:n
      j, dj = jacobi_and_derivative(ξ,i-3,1,1)
      _, d2j = jacobi_and_derivative(ξ,i-4,2,2)
      @inbounds v[d,i] = -sqrt(2*i-3)*(2*dv1*dv2*j+2*(dv1*v2+v1*dv2)*(2/(b[d]-a[d]))*dj+v1*v2*d2j*2*i*((b[d]-a[d])^2))/(i-2)
    end
  end
end

function _evaluate_nd_mc0!(
  v::AbstractVector{V},
  x,
  a::Vector{Point{D,T}},
  b::Vector{Point{D,T}},
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  o = one(T)
  k = 1
  l = length(terms)

  for (i,ci) in enumerate(terms)

    for d in 1:dim
      _evaluate_1d_mc0!(c,x,a[i],b[i],orders[d],d)
    end

    s = o
    for d in 1:dim
      @inbounds s *= c[d,ci[d]]
    end

    k = _set_value_mc0!(v,s,k,l)

  end

end

@inline function _set_value_mc0!(v::AbstractVector{V},s::T,k,l) where {V,T}
  m = zero(Mutable(V))
  z = zero(T)
  js = eachindex(m)
  for j in js
    for i in js
      @inbounds m[i] = z
    end
    @inbounds m[j] = s
    i = k+l*(j-1)
    @inbounds v[i] = m
  end
  k+1
end

@inline function _set_value_mc0!(v::AbstractVector{<:Real},s,k,l)
  @inbounds v[k] = s
  k+1
end

function _gradient_nd_mc0!(
  v::AbstractVector{G},
  x,
  a::Vector{Point{D,T}},
  b::Vector{Point{D,T}},
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  z = zero(Mutable(VectorValue{D,T}))
  o = one(T)
  k = 1
  l = length(terms)

  for (i,ci) in enumerate(terms)

    for d in 1:dim
      _evaluate_1d_mc0!(c,x,a[i],b[i],orders[d],d)
      _gradient_1d_mc0!(g,x,a[i],b[i],orders[d],d)
    end

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

    k = _set_gradient_mc0!(v,s,k,l,V)

  end

end

@inline function _set_gradient_mc0!(
  v::AbstractVector{G},s,k,l,::Type{<:Real}) where G

  @inbounds v[k] = s
  k+1
end

@inline function _set_gradient_mc0!(
  v::AbstractVector{G},s,k,l,::Type{V}) where {V,G}

  T = eltype(s)
  m = zero(Mutable(G))
  w = zero(V)
  z = zero(T)
  for (ij,j) in enumerate(CartesianIndices(w))
    for i in CartesianIndices(m)
      @inbounds m[i] = z
    end
    for i in CartesianIndices(s)
      @inbounds m[i,j] = s[i]
    end
    i = k+l*(ij-1)
    @inbounds v[i] = m
  end
  k+1
end

function _hessian_nd_mc0!(
  v::AbstractVector{G},
  x,
  a::Vector{Point{D,T}},
  b::Vector{Point{D,T}},
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  z = zero(Mutable(TensorValue{D,D,T}))
  o = one(T)
  k = 1
  l = length(terms)

  for (i,ci) in enumerate(terms)

    for d in 1:dim
      _evaluate_1d_mc0!(c,x,a[i],b[i],orders[d],d)
      _gradient_1d_mc0!(g,x,a[i],b[i],orders[d],d)
      _hessian_1d_mc0!(h,x,a[i],b[i],orders[d],d)
    end

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

    k = _set_gradient_mc0!(v,s,k,l,V)

  end

end
