"""
    ModalC0 <: Polynomial

Type representing ModalC0 polynomials, c.f. [ModalC0 polynomials](@ref) section.

Reference: Eq. (17) in https://doi.org/10.1016/j.camwa.2022.09.027
"""
struct ModalC0 <: Polynomial end

"""
    ModalC0Basis{D,V,T,K} <: PolynomialBasis{D,V,K,ModalC0}

Tensor product basis of generalised modal C0 1D basis from section 5.2 in
https://doi.org/10.1016/j.camwa.2022.09.027.
See also [ModalC0 polynomials](@ref) section of the documentation.
"""
struct ModalC0Basis{D,V,T,K} <: PolynomialBasis{D,V,K,ModalC0}
  orders::NTuple{D,Int}
  terms::Vector{CartesianIndex{D}}
  a::Vector{Point{D,T}}
  b::Vector{Point{D,T}}

  function ModalC0Basis{D}(
    ::Type{V},
    orders::NTuple{D,Int},
    terms::Vector{CartesianIndex{D}},
    a::Vector{Point{D,T}},
    b::Vector{Point{D,T}}) where {D,V,T}

    @check T == eltype(V) "Point and polynomial values should have the same scalar body"
    K = maximum(orders, init=0)

    new{D,V,T,K}(orders,terms,a,b)
  end
end

function ModalC0Basis{D}(
  ::Type{V},
  orders::NTuple{D,Int},
  a::Vector{Point{D,T}},
  b::Vector{Point{D,T}};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,V,T}

  terms = _define_terms_mc0(filter, sort!, orders)
  ModalC0Basis{D}(V,orders,terms,a,b)
end

function ModalC0Basis{D}(
  ::Type{V},
  orders::NTuple{D,Int},
  sa::Point{D,T},
  sb::Point{D,T};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,V,T}

  terms = _define_terms_mc0(filter, sort!, orders)
  a = fill(sa,length(terms))
  b = fill(sb,length(terms))
  ModalC0Basis{D}(V,orders,terms,a,b)
end

function ModalC0Basis{D}(
  ::Type{V},
  orders::NTuple{D,Int};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,V}

  T = eltype(V)
  sa = Point{D,T}(tfill(zero(T),Val{D}()))
  sb = Point{D,T}(tfill( one(T),Val{D}()))
  ModalC0Basis{D}(V,orders,sa,sb; filter=filter, sort! =sort!)
end

function ModalC0Basis{D}(
  ::Type{V},
  order::Int,
  a::Vector{Point{D,T}},
  b::Vector{Point{D,T}};
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,V,T}

  orders = tfill(order,Val{D}())
  ModalC0Basis{D}(V,orders,a,b; filter=filter, sort! =sort!)
end

function ModalC0Basis{D}(
  ::Type{V},
  order::Int;
  filter::Function=_q_filter,
  sort!::Function=_sort_by_nfaces!) where {D,V}

  orders = tfill(order,Val{D}())
  ModalC0Basis{D}(V,orders; filter=filter, sort! =sort!)
end

# API

@inline Base.size(a::ModalC0Basis{D,V}) where {D,V} = (length(a.terms)*num_indep_components(V),)

function get_orders(b::ModalC0Basis)
  b.orders
end

# Helpers

_s_filter_mc0(e,o) = ( sum( [ i for i in e if i>1 ] ) <= o )

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
  maxorder = maximum(orders)
  term_to_is_fterm = lazy_map(t->filter(Int[Tuple(t-to)...],maxorder),terms)
  findall(term_to_is_fterm)
end

function _define_terms_mc0(filter,sort!,orders)
  terms = _define_terms(_q_filter,orders)
  sort!(terms,orders)
  mask = _compute_filter_mask(terms,filter,orders)
  collect(lazy_map(Reindex(terms),mask))
end


#################################
# nD evaluations implementation #
#################################

function _evaluate_nd!(
  basis::ModalC0Basis{D,V,T,K}, x,
  r::AbstractMatrix{V}, i,
  c::AbstractMatrix{T}) where {D,V,T,K}

  terms  = basis.terms
  orders = basis.orders
  a = basis.a
  b = basis.b

  k = 1
  l = length(terms)

  for (n,ci) in enumerate(terms)

    for d in 1:D
      _evaluate_1d_mc0!(c,x,a[n],b[n],orders[d],d)
    end

    s = one(T)
    for d in 1:D
      @inbounds s *= c[d,ci[d]]
    end

    k = _set_value_mc0!(r,i,s,k,l)
  end
end

@inline function _set_value_mc0!(r::AbstractMatrix{V},i,s::T,k,l) where {V,T}
  ncomp = num_indep_components(V)
  z = zero(T)
  for j in 1:ncomp
    m = k+l*(j-1)
    @inbounds r[i,m] = ntuple(p -> ifelse(p == j, s, z),Val(ncomp))
  end
  k+1
end

@inline function _set_value_mc0!(r::AbstractMatrix{<:Real},i,s,k,l)
  @inbounds r[i,k] = s
  k+1
end

function _gradient_nd!(
  basis::ModalC0Basis{D,V,T,K}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}) where {D,V,T,K,G}

  terms  = basis.terms
  orders = basis.orders
  a = basis.a
  b = basis.b

  k = 1
  l = length(terms)

  for (n,ci) in enumerate(terms)

    for d in 1:D
      _evaluate_1d_mc0!(c,x,a[n],b[n],orders[d],d)
      _gradient_1d_mc0!(g,x,a[n],b[n],orders[d],d)
    end

    for j in eachindex(s)
      @inbounds s[j] = one(T)
    end
    for q in 1:D
      for d in 1:D
        if d != q
          @inbounds s[q] *= c[d,ci[d]]
        else
          @inbounds s[q] *= g[d,ci[d]]
        end
      end
    end

    k = _set_derivative_mc0!(r,i,s,k,l,V)
  end
end

@inline function _set_derivative_mc0!(
  r::AbstractMatrix{G},i,s,k,l,::Type{<:Real}) where G

  @inbounds r[i,k] = s
  k+1
end

# Indexing and m definition should be fixed if G contains symmetries, that is
# if the code is  optimized for symmetric tensor V valued FESpaces
# (if gradient_type(V) returned a symmetric higher order tensor type G)
@inline @generated function _set_derivative_mc0!(
  r::AbstractMatrix{G},i1,s,k,l,::Type{V}) where {V,G}
  # Git blame me for readable non-generated version
  @notimplementedif num_indep_components(V) != num_components(V) "Not implemented for symmetric Jacobian or Hessian"

  m = Array{String}(undef, size(G))
  N_val_dims = length(size(V))
  s_size = size(G)[1:end-N_val_dims]

  body = "T = eltype(s); z = zero(T);"
  for ci in CartesianIndices(s_size)
    id = join(Tuple(ci))
    body *= "@inbounds s$id = s[$ci];"
  end

  V_size = size(V)
  for (ij,j) in enumerate(CartesianIndices(V_size))
    for i in CartesianIndices(m)
      m[i] = "z"
    end
    for ci in CartesianIndices(s_size)
      id = join(Tuple(ci))
      m[ci,j] = "s$id"
    end
    body *= "i = k + l*($ij-1);"
    body *= "@inbounds r[i1,i] = ($(join(tuple(m...), ", ")));"
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k+1))
end

function _hessian_nd!(
  basis::ModalC0Basis{D,V,T,K}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  s::MMatrix{D,D,T}) where {D,V,T,K,G}

  terms  = basis.terms
  orders = basis.orders
  a = basis.a
  b = basis.b

  k = 1
  l = length(terms)

  for (n,ci) in enumerate(terms)

    for d in 1:D
      _evaluate_1d_mc0!(c,x,a[n],b[n],orders[d],d)
      _gradient_1d_mc0!(g,x,a[n],b[n],orders[d],d)
      _hessian_1d_mc0!(h,x,a[n],b[n],orders[d],d)
    end

    for j in eachindex(s)
      @inbounds s[j] = one(T)
    end
    for r in 1:D
      for q in 1:D
        for d in 1:D
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

    k = _set_derivative_mc0!(r,i,s,k,l,V)
  end
end


#################################
# 1D evaluations implementation #
#################################

# Reference: equation (17) in
#
# Badia, S.; Neiva, E. & Verdugo, F.; (2022);
# Robust high-order unfitted finite elements by interpolation-based discrete extension,
# Computers & Mathematics with Applications,
# https://doi.org/10.1016/j.camwa.2022.09.027
function _evaluate_1d_mc0!(c::AbstractMatrix{T},x,a,b,order,d) where T
  @assert order > 0
  n = order + 1
  o = one(T)
  @inbounds c[d,1] = o - x[d]
  @inbounds c[d,2] = x[d]
  if n > 2
    ξ = ( 2*x[d] - ( a[d] + b[d] ) ) / ( b[d] - a[d] )
    for i in 3:n
      @inbounds c[d,i] = -sqrt(2*i-3)*c[d,1]*c[d,2]*jacobi(ξ,i-3,1,1)/(i-2)
    end
  end
end

function _gradient_1d_mc0!(g::AbstractMatrix{T},x,a,b,order,d) where T
  @assert order > 0
  n = order + 1
  z = one(T)
  @inbounds g[d,1] = -z
  @inbounds g[d,2] = z
  if n > 2
    ξ = ( 2*x[d] - ( a[d] + b[d] ) ) / ( b[d] - a[d] )
    v1 = z - x[d]
    v2 = x[d]
    for i in 3:n
      j, dj = jacobi_and_derivative(ξ,i-3,1,1)
      @inbounds g[d,i] = -sqrt(2*i-3)*(g[d,1]*v2*j+v1*g[d,2]*j+v1*v2*(2/(b[d]-a[d]))*dj)/(i-2)
    end
  end
end

function _hessian_1d_mc0!(h::AbstractMatrix{T},x,a,b,order,d) where T
  @assert order > 0
  n = order + 1
  y = zero(T)
  z = one(T)
  @inbounds h[d,1] = y
  @inbounds h[d,2] = y
  if n > 2
    ξ = ( 2*x[d] - ( a[d] + b[d] ) ) / ( b[d] - a[d] )
    v1 = z - x[d]
    v2 = x[d]
    dv1 = -z
    dv2 = z
    for i in 3:n
      j, dj = jacobi_and_derivative(ξ,i-3,1,1)
      _, d2j = jacobi_and_derivative(ξ,i-4,2,2)
      @inbounds h[d,i] = -sqrt(2*i-3)*(2*dv1*dv2*j+2*(dv1*v2+v1*dv2)*(2/(b[d]-a[d]))*dj+v1*v2*d2j*2*i*((b[d]-a[d])^2))/(i-2)
    end
  end
end

