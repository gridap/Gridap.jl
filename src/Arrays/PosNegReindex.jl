
# """
# """
# struct LocalToGlobalPosNegArray{T,M,N,L,V,W} <: AbstractArray{Array{T,M},N}
#   lid_to_gid::L
#   gid_to_val_pos::V
#   gid_to_val_neg::W

#   @doc """
#   """
#   function LocalToGlobalPosNegArray(
#     lid_to_gid::AbstractArray{<:AbstractArray},
#     gid_to_val_pos::AbstractArray,
#     gid_to_val_neg::AbstractArray)

#     L = typeof(lid_to_gid)
#     V = typeof(gid_to_val_pos)
#     W = typeof(gid_to_val_neg)
#     T = eltype(V)
#     M = ndims(eltype(L))
#     N = ndims(L)
#     new{T,M,N,L,V,W}(lid_to_gid,gid_to_val_pos,gid_to_val_neg)
#   end
# end

# size(a::LocalToGlobalPosNegArray) = size(a.lid_to_gid)

# function IndexStyle(::Type{<:LocalToGlobalPosNegArray{T,M,N,L}}) where {T,M,N,L}
#   IndexStyle(L)
# end

# function array_cache(a::LocalToGlobalPosNegArray)
#   gids = testitem(a.lid_to_gid)
#   T = eltype(a.gid_to_val_pos)
#   r = zeros(T,size(gids))
#   c = CachedArray(r)
#   cl = array_cache(a.lid_to_gid)
#   (cl,c)
# end

# function getindex!(cache,a::LocalToGlobalPosNegArray,i::Integer...)
#   (cl,c) = cache
#   gids = getindex!(cl,a.lid_to_gid,i...)
#   setsize!(c,size(gids))
#   r = c.array
#   for (i,gid) in enumerate(gids)
#     if gid > 0
#       r[i] = a.gid_to_val_pos[gid]
#     elseif gid < 0
#       r[i] = a.gid_to_val_neg[-gid]
#     else
#       @unreachable "Only positive or negative indices allowed, not zero."
#     end
#   end
#   r
# end

# function getindex(a::LocalToGlobalPosNegArray,i::Integer...)
#   cache = array_cache(a)
#   getindex!(cache,a,i...)
# end

# This Map has non-trivial domain, thus we need the define testargs
"""
    PosNegReindex(values_pos,values_neg)
"""
struct PosNegReindex{A,B} <: Map
  values_pos::A
  values_neg::B
end

function testargs(k::PosNegReindex,i::Integer)
  @check length(k.values_pos) !=0 || length(k.values_neg) != 0 "This map has empty domain"
  @check eltype(k.values_pos) == eltype(k.values_neg) "This map is type-instable"
  length(k.values_pos) !=0 ? (one(i),) : (-one(i))
end

function return_value(k::PosNegReindex,i::Integer)
  if length(k.values_pos)==0 && length(k.values_neg)==0
    @check eltype(k.values_pos) == eltype(k.values_neg) "This map is type-instable"
    testitem(k.values_pos)
  else
    evaluate(k,testargs(k,i)...)
  end
end

function return_cache(k::PosNegReindex,i::Integer)
  c_p = array_cache(k.values_pos)
  c_n = array_cache(k.values_neg)
  c_p, c_n
end

@inline function evaluate!(cache,k::PosNegReindex,i::Integer)
  c_p, c_n = cache
  i>0 ? getindex!(c_p,k.values_pos,i) : getindex!(c_n,k.values_neg,-i)
end


#@inline return_type(k::PosNegReindex,x...) = typeof(first(k.values))
#
#@inline return_type(k::PosNegReindex,x::AbstractArray...) = typeof(testitem(k,x...))
#
#@inline return_cache(k::PosNegReindex,f) = array_cache(k.values)
#
#@inline function return_cache(k::PosNegReindex,a::AbstractArray)
#  gids = a
#  values_pos = k.values_pos
#  values_neg = k.values_neg
#  T = eltype(values_pos)
#  S = eltype(values_neg)
#  @assert T == S "The types of both value arrays must be the same"
#  r = zeros(T,size(gids))
#  c = CachedArray(r)
#end
#
#@inline function evaluate!(cache,k::PosNegReindex,i)
#  i > 0 ? getindex!(cache,k.values_pos,i) : getindex!(cache,k.values_neg,i)
#end
#
#@inline function evaluate!(cache,k::PosNegReindex,gids::AbstractArray)
#  setsize!(cache,size(gids))
#  r = cache.array
#  @inbounds for (i,gid) in enumerate(gids)
#    if gid > 0
#      r[i] = k.values_pos[gid]
#    elseif gid < 0
#      r[i] = k.values_neg[-gid]
#    else
#      @unreachable "Only positive or negative indices allowed, not zero."
#    end
#  end
#  r
#end

#function posneg_reindex(i_to_v_pos::AbstractArray, i_to_v_neg::AbstractArray, j_to_i::AbstractArray)
#  lazy_map(PosNegReindex(i_to_v_pos,i_to_v_neg),j_to_i)
#end
#
#function posneg_reindex(i_to_v_pos::AbstractArray, i_to_v_neg::AbstractArray, j_to_i::AbstractArray{<:AbstractArray})
#  lazy_map(PosNegReindex(i_to_v_pos,i_to_v_neg),j_to_i)
#end

# @propagate_inbounds function Base.setindex!(a::LazyArray{<:Fill{<:PosNegReindex}},v,j::Integer)
#   k = a.g.value
#   i_to_v = k.values
#   j_to_i, = a.f
#   i = j_to_i[j]
#   i_to_v[i]=v
# end

#@inline function testitem(k::PosNegReindex,gids)
#  gids == 0 ? testitem(k,1) :  evaluate(k,gids)#testvalue(eltype(k.values))
#end
