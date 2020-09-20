# Apply an array of Mapping to other arrays

function apply(g::AbstractArray,f::AbstractArray...)
  MappedArray(g,f...)
end

function apply(::Type{T},g::AbstractArray,f::AbstractArray...) where T
  MappedArray(T,g,f...)
end

function apply_mapping(k,f::AbstractArray...)
    s = common_size(f...)
    apply(Fill(k, s...), f...)
end

function apply_mapping(::Type{T},k,f::AbstractArray...) where T
  s = common_size(f...)
  apply(T,Fill(k, s...), f...)
end

apply(T::Type,k::Mapping,f::AbstractArray...) = apply_mapping(T,k,f...)

apply(T::Type,k::Function,f::AbstractArray...) = apply_mapping(T,k,f...)

apply(k::Mapping,f::AbstractArray...) = apply_mapping(k,f...)

apply(k::Function,f::AbstractArray...) = apply_mapping(k,f...)

struct MappedArray{G,T,N,F} <: AbstractArray{T,N}
  g::G
  f::F
  function MappedArray(::Type{T}, g::AbstractArray, f::AbstractArray...) where T
    G = typeof(g)
    F = typeof(f)
    new{G,T,ndims(first(f)),F}(g, f)
  end
end

function MappedArray(g::AbstractArray, f::AbstractArray...)
  gi = testitem(g) # Assumes that all mappings return the same type
  fi = testitems(f...)
  T = typeof(testitem(gi, fi...))
  MappedArray(T, g, f...)
end

IndexStyle(a::MappedArray) = IndexStyle(a.g)

# @fverdugo : TODO a lot of splating i... which is bad for CartesianIndices

function array_cache(a::MappedArray,i...)
  @notimplementedif ! all(map(isconcretetype, map(eltype, a.f)))
  if ! (eltype(a.g) <: Function)
    @notimplementedif ! isconcretetype(eltype(a.g))
  end
  gi = testitem(a.g)
  fi = testitems(a.f...)
  cg = return_cache(a.g,i...)
  cf = return_caches(a.f,i...)
  cgi = return_cache(gi, fi...)
  cg, cgi, cf
end

function array_cache(a::MappedArray)
@notimplementedif ! all(map(isconcretetype, map(eltype, a.f)))
if ! (eltype(a.g) <: Function)
  @notimplementedif ! isconcretetype(eltype(a.g))
end
gi = testitem(a.g)
fi = testitems(a.f...)
cg = array_cache(a.g)
cf = map(array_cache,a.f)
cgi = return_cache(gi, fi...)
cg, cgi, cf
end

@inline function getindex!(cache, a::MappedArray, i...)
  cg, cgi, cf = cache
  gi = getindex!(cg, a.g, i...)
  fi = getitems!(cf, a.f, i...)
  vi = evaluate!(cgi, gi, fi...)
  vi
end

function Base.getindex(a::MappedArray, i...)
  ca = array_cache(a)
  getindex!(ca, a, i...)
end

function Base.size(a::MappedArray)
  size(first(a.f))
end

# Particular implementations for Fill

function apply_mapping(f::Fill, a::Fill...)
  ai = _getvalues(a...)
  r = evaluate(f.value, ai...)
  s = common_size(f, a...)
  Fill(r, s)
end

function apply_mapping(::Type{T}, f::Fill, a::Fill...) where T
  apply_mapping(f, a...)
end

function _getvalues(a::Fill, b::Fill...)
  ai = a.value
  bi = _getvalues(b...)
  (ai, bi...)
end

function _getvalues(a::Fill)
  ai = a.value
  (ai,)
end

# Operator

function apply(op::MappingOperator,x...)
  apply(Fill(op,length(first(x))),x...)
end

function apply(
  ::typeof(evaluate),
  a::MappedArray{<:Fill{<:MappingOperator}},
  x::AbstractArray)

  fx = map( fi->apply(evaluate,fi,x), a.f)
  op = a.g.value.op
  apply( (args...) -> op(args...), fx...)
  # apply( (args...) -> broadcast(op,args...), fx... )  # TODO use BroadcastMapping
end

# function return_mapping_array_cache(a::AbstractArray, x::AbstractArray)
#   ca = array_cache(a)
#   fi = testitem(a)
#   xi = testitem(x)
#   cfi = return_cache(fi, xi)
#   cx = array_cache(x)
#   (ca, cfi, cx)
# end

function test_mapped_array(
  a::AbstractArray,
  x::AbstractArray,
  v::AbstractArray,
  cmp::Function=(==);
  grad=nothing)

  ax = apply_mapping(a, x)
  test_array(ax, v, cmp)

  ca, cfi, cx = array_cache(a, x)
  t = true
  for i in 1:length(a)
    fi = getindex!(ca, a, i)
    xi = getindex!(cx, x, i)
    fxi = evaluate!(cfi, fi, xi)
    vi = v[i]
    ti = cmp(fxi, vi)
    t = t && ti
  end
  @test t

  if grad != nothing
    g = apply_mapping(gradient, a)
    test_mapped_array(g, x, grad, cmp)
  end

end

function common_size(a::AbstractArray...)
  a1, = a
  c = all([size(a1) == size(ai) for ai in a])
  if !c
    error("Array sizes $(map(size,a)) are not compatible.")
  end
  l = size(a1)
  (l,)
end
