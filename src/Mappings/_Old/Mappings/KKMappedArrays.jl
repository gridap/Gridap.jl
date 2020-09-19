# @santiagobadia : Change names later on apply...
apply_mapping(f::Function,a...) = _apply_broadcast(f,a...)

apply_mapping(f::Mapping,a...) = _apply_broadcast(f,a...)

function _apply_broadcast(f,a::AbstractArray...)
    s = common_size(a...)
    apply_mapping(Fill(f, s...), a...)
end

function apply_mapping(::Type{T}, f, a::AbstractArray...) where T
    s = common_size(a...)
    apply_mapping(T, Fill(f, s...), a...)
end

function apply_mapping(f::AbstractArray, a::AbstractArray...)
    MappedArray(f, a...)
end

function apply_mapping(::Type{T}, f::AbstractArray, a::AbstractArray...) where T
    MappedArray(T, f, a...)
end

struct MappedArray{G,T,N,F} <: AbstractArray{T,N}
    g::G
    f::F
    function MappedArray(::Type{T}, g::AbstractArray, f::AbstractArray...) where T
        G = typeof(g)
        F = typeof(f)
        f1, = f
        new{G,T,ndims(f1),F}(g, f)
    end
end

function MappedArray(g::AbstractArray, f::AbstractArray...)
    gi = testitem(g) # Assumes that all mappings return the same type
    fi = testitems(f...)
    T = typeof(testitem(gi, fi...))
    MappedArray(T, g, f...)
end

IndexStyle(a::MappedArray) = IndexStyle(a.g)
# function uses_hash(::Type{<:MappedArray})
#     Val(true)
# end

# function array_cache(hash::Dict, a::MappedArray)
#     id = objectid(a)
#     _cache = _array_cache(hash, a)
#     if haskey(hash, id)
#         cache = _get_stored_cache(_cache, hash, id)
#     else
#         cache = _cache
#         hash[id] = cache
#     end
#     cache
# end

# # TODO
# @noinline function _get_stored_cache(cache::T, hash, id) where T
#     hash[id]
# end

# function _array_cache(hash, a::MappedArray)
#     cg = array_cache(hash, a.g)
#     gi = testitem(a.g)
#     fi = testitems(a.f...)
#     @notimplementedif ! all(map(isconcretetype, map(eltype, a.f)))
#     if ! (eltype(a.g) <: Function)
#         @notimplementedif ! isconcretetype(eltype(a.g))
#     end
#     cf = array_caches(hash, a.f...)
#     cgi = return_cache(gi, fi...)
#     ai = testitem!(cgi, gi, fi...)
#     i = -testitem(eachindex(a))
#     e = Evaluation((i,), ai)
#     c = (cg, cgi, cf)
#     (c, e)
# end

function testitem(a::MappedArray)
    # cg = array_cache(a.g)
    gi = testitem(a.g)
    fi = testitems(a.f...)
    testitem(gi, fi...)
end

# function getindex!(cache, a::MappedArray, i::Integer...)
#     li = LinearIndices(a)
#     getindex!(cache, a, li[i...])
# end

# function getindex!(cache, a::MappedArray, i::Integer)
#     _cached_getindex!(cache, a, (i,))
# end

# function getindex!(cache, a::MappedArray, i::CartesianIndex)
#     _cached_getindex!(cache, a, Tuple(i))
# end

# function _cached_getindex!(cache, a::MappedArray, i::Tuple)
#     c, e = cache
#     v = e.fx
#     if e.x != i
#         v = _getindex!(c, a, i...)
#         e.x = i
#         e.fx = v
#     end
#     v
# end

function _getindex!(cache, a::MappedArray, i...)
    cg, cgi, cf = cache
    gi = getindex!(cg, a.g, i...)
    fi = getitems!(cf, a.f, i...)
    vi = evaluate!(cgi, gi, fi...)
    vi
end

function Base.getindex(a::MappedArray, i...)
    ca = return_cache(a)
    getindex!(ca, a, i...)
end

function Base.IndexStyle(
  ::Type{MappedArray{T,N,F,G}}) where {T,N,F,G}
    common_index_style(F)
end

function Base.size(a::MappedArray)
    f, = a.f
    size(f)
end

# function common_size(a::AbstractArray...)
#     a1, = a
#     c = all([length(a1) == length(ai) for ai in a])
#     if !c
#         error("Array sizes $(map(size, a)) are not compatible.")
#     end
#     l = length(a1)
#     (l,)
# end

# # TODO not sure what to do with shape and index-style
# function common_index_style(::Type{<:Tuple})
#     IndexLinear()
# end

# mutable struct Evaluation{X,F} <: GridapType
#     x::X
#     fx::F
#     function Evaluation(x::X, fx::F) where {X,F}
#         new{X,F}(x, fx)
#     end
# end

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

# function Base.sum(a::MappedArray)
#     cache = array_cache(a)
#     _sum_array_cache(cache, a)
# end

# function _sum_array_cache(cache, a)
#     r = zero(eltype(a))
#     for i in eachindex(a)
#         ai = getindex!(cache, a, i)
#         r += ai
#     end
#     r
# end

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

  ca, cfi, cx = return_mapping_array_cache(a, x)
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
