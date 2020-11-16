
"""
"""
function lazy_append(a::AbstractArray,b::AbstractArray)
  AppendedArray(a,b)
end

"""
"""
function lazy_split(f::AbstractArray,n::Integer)
  _lazy_split(f,n)
end

function lazy_split(f::CompressedArray,n::Integer)
  _a , _b =_lazy_split(f,n)
  a = _compact_values_ptrs(_a)
  b = _compact_values_ptrs(_b)
  a,b
end

function _lazy_split(f,n)
  l = length(f)
  @assert n <= l
  ids_a = collect(1:n)
  ids_b = collect((n+1):l)
  a = lazy_map(Reindex(f),ids_a)
  b = lazy_map(Reindex(f),ids_b)
  a,b
end

function _compact_values_ptrs(a)
  v = a.values
  p = a.ptrs
  touched = fill(false,length(v))
  for i in p
    touched[i] = true
  end
  if all(touched) || length(p) == 0
    return a
  else
    j_to_i = findall(touched)
    i_to_j = zeros(Int,length(v))
    i_to_j[j_to_i] = 1:length(j_to_i)
    _v = v[j_to_i]
    _p = i_to_j[p]
    return CompressedArray(_v, _p)
  end
end

struct AppendedArray{T,A,B} <: AbstractVector{T}
  a::A
  b::B
  function AppendedArray(a::AbstractArray,b::AbstractArray)
    A = typeof(a)
    B = typeof(b)
    ai = testitem(a)
    bi = testitem(b)
    T = eltype(collect((ai,bi)))
    new{T,A,B}(a,b)
  end
end

Base.IndexStyle(::Type{<:AppendedArray}) = IndexLinear()

function Base.getindex(v::AppendedArray,i::Integer)
  l = length(v.a)
  if i > l
    v.b[i-l]
  else
    v.a[i]
  end
end

Base.size(v::AppendedArray) = (length(v.a)+length(v.b),)

function array_cache(v::AppendedArray)
  ca = array_cache(v.a)
  cb = array_cache(v.b)
  (ca,cb)
end

function getindex!(cache,v::AppendedArray,i::Integer)
  ca, cb = cache
  l = length(v.a)
  if i > l
    getindex!(cb,v.b,i-l)
  else
    getindex!(ca,v.a,i)
  end
end

#function lazy_map(f,a::AppendedArray...)
#  la = map(ai->length(ai.a),a)
#  lb = map(ai->length(ai.b),a)
#  if all(la .== first(la)) && all(lb .== first(lb))
#    c_a = lazy_map(f,map(ai->ai.a,a)...)
#    c_b = lazy_map(f,map(ai->ai.b,a)...)
#    lazy_append(c_a,c_b)
#  else
#    s = _common_size(a...)
#    LazyArray(Fill(f,s...),a...)
#  end
#end
#
#function lazy_map(f,::Type{T},a::AppendedArray...) where T
#  la = map(ai->length(ai.a),a)
#  lb = map(ai->length(ai.b),a)
#  if all(la .== first(la)) && all(lb .== first(lb))
#    c_a = lazy_map(f,T,map(ai->ai.a,a)...)
#    c_b = lazy_map(f,T,map(ai->ai.b,a)...)
#    lazy_append(c_a,c_b)
#  else
#    s = _common_size(a...)
#    LazyArray(Fill(f,s...),T,a...)
#  end
#end

function lazy_map(::typeof(evaluate),::Type{T},b::Fill,a::AppendedArray...) where T
  f = b.value
  la = map(ai->length(ai.a),a)
  lb = map(ai->length(ai.b),a)
  if all(la .== first(la)) && all(lb .== first(lb))
    c_a = lazy_map(f,map(ai->ai.a,a)...)
    c_b = lazy_map(f,map(ai->ai.b,a)...)
    lazy_append(c_a,c_b)
  else
    s = _common_size(a...)
    LazyArray(T,Fill(f,s...),a...)
  end
end

function lazy_map(k::typeof(evaluate),::Type{T},a::AppendedArray...) where T
  @check all(map(i->length(i)==length(a[1]),a))
  la = map(ai->length(ai.a),a)
  if all(la .== first(la))
    c_a = lazy_map(k,map(ai->ai.a,a)...)
    c_b = lazy_map(k,map(ai->ai.b,a)...)
    lazy_append(c_a,c_b)
  else
    LazyArray(T,a...)
  end
end

# Very important optimization to compute error norms efficiently e.g. in EmbeddedFEM
function Base.sum(a::AppendedArray)
  sum(a.a) + sum(a.b)
end

