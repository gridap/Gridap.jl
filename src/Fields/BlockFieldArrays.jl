
function similar_range(r::Base.OneTo,n::Integer)
  Base.OneTo(Int(n))
end

function similar_range(r::BlockedUnitRange,n::Integer)
  blockedrange([n])
end

function similar_range(r::MultiLevelBlockedUnitRange,n::Integer)
  r = similar_range(first(r.local_ranges),n)
  append_ranges([r])
end

struct BlockFieldArrayCoo{T,N,A,X} <: AbstractBlockArray{T,N}
  axes::X
  blockids::Vector{NTuple{N,Int}}
  blocks::Vector{A}
  ptrs::Array{Int,N}
  function BlockFieldArrayCoo(
    axes::NTuple{N},
    blockids::Vector{NTuple{N,Int}},
    blocks::Vector{A},
    ptrs::Array{Int,N}) where {T,N,A<:AbstractArray{T,N}}

    @assert length(blockids) == length(blocks)
    @check length(unique(blockids)) == length(blockids) "We cannot built a BlockFieldArrayCoo from repeated blocks"
    @check Arrays._valid_block_sizes(axes,ptrs,blocks) "The given blocks do not match with the given axes"

    X = typeof(axes)
    new{T,N,A,X}(axes,blockids,blocks,ptrs)
  end
end

# Minimal constructor

function BlockFieldArrayCoo(
  axes::NTuple{N},
  blockids::Vector{NTuple{N,Int}},
  blocks::Vector{A}) where {T,N,A<:AbstractArray{T,N}}

  ptrs = Arrays._compute_ptrs(blockids,axes)
  BlockFieldArrayCoo(axes,blockids,blocks,ptrs)
end

# Minimal constructor (for lazy_map)

function BlockFieldArrayCoo(
  axes::NTuple{N},
  blockids::Vector{NTuple{N,Int}},
  blocks::A...) where {T,N,A<:AbstractArray{T,N}}

  BlockFieldArrayCoo(axes,blockids,collect(blocks))
end

function return_cache(
  ::Type{<:BlockFieldArrayCoo},
  axes::NTuple{N},
  blockids::Vector{NTuple{N,Int}},
  blocks::A...) where {T,N,A<:AbstractArray{T,N}}

  r = BlockFieldArrayCoo(axes,blockids,blocks...)
  CachedArray(r)
end

@inline function evaluate!(
  cache,
  ::Type{<:BlockFieldArrayCoo},
  axes::NTuple{N},
  blockids::Vector{NTuple{N,Int}},
  blocks::A...) where {T,N,A<:AbstractArray{T,N}}
                           
  setaxes!(cache,axes)
  r = cache.array
  copyto!(r.blocks,blocks)
  r
end

# Specific API

function is_zero_block(a::BlockFieldArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  p = a.ptrs[i...]
  @check p != 0
  p < 0
end

# AbstractBlockArray

@inline function BlockArrays.getblock(a::BlockFieldArrayCoo{T,N}, block::Vararg{Integer, N}) where {T,N}
  #@boundscheck BlockArrays.blockcheckbounds(a, block...)
  p = a.ptrs[block...]
  if p>0
    a.blocks[p]
  else
    @notimplemented "Cannot get a zero block from a BlockFieldArrayCoo"
  end
end

function BlockArrays.getblock!(c,a::BlockFieldArrayCoo{T,N}, block::Vararg{Integer, N}) where {T,N}
  b = a[Block(block...)]
  copy!(c,b)
  c
end

# AbstractArray

Base.size(a::BlockFieldArrayCoo) = map(length,Base.axes(a))
Base.axes(a::BlockFieldArrayCoo) = a.axes
Base.IndexStyle(::Type{<:BlockFieldArrayCoo}) = IndexCartesian()

function Base.getindex(a::BlockFieldArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

# Evaluation

function return_cache(a::BlockFieldArrayCoo,x::Point)
  fc = map(i->return_cache(i,x),a.blocks)
  fx = map(i->return_value(i,x),a.blocks)
  cr = return_cache(BlockArrayCoo,a.axes,a.blockids,fx...)
  (fc,fx,cr)
end

@inline function evaluate!(cache,a::BlockFieldArrayCoo,x::Point)
  fc,fx, cr = cache
  #for i in 1:length(a.blocks)
  #  fx[i] = evaluate!(cr[i],a.blocks[i],x)
  #end
  #evaluate!(cr,BlockArrayCoo,a.axes,a.blockids,fx...)
  cr.array
end

## The following masted types are needed to achieve type-stability
## in order to use BlockVectorCoo with arrays of fields
#
#struct MaskedField{F} <: Field
#  field::F
#  mask::Bool
#end
#
#function return_cache(z::MaskedField,x::Point)
#  return_cache(z.field,x)
#end
#
#@inline function evaluate!(cache,z::MaskedField,x::Point) 
#  fx = evaluate!(cache,z.field,x)
#  z.mask ? zero(fx) : fx
#end
#
#testvalue(::Type{MaskedField{F}}) where F = MaskedField(testvalue(F),false)
#
#function return_cache(z::MaskedField,x::AbstractArray{<:Point})
#  return_cache(z.field,x)
#end
#
#function evaluate!(c,z::MaskedField,x::AbstractArray{<:Point})
#  fx = evaluate!(c,z.field,x)
#  if z.mask
#    fill!(fx,zero(eltype(fx)))
#  end
#  fx
#end
#
#@inline gradient(z::MaskedField) = MaskedField(gradient(z.field),z.mask)
#
#struct MaskedFieldArray{T,N,A,X} <: AbstractArray{T,N}
#  axes::X
#  field_array::A
#  mask::Bool
#  function MaskedFieldArray(
#    _axes::NTuple{N},
#    field_array::AbstractArray{F,N},
#    mask::Bool) where {F<:Field,N}
#
#    @check mask || blocks_equal(_axes,axes(field_array)) "Incompatible axes and field_array"
#
#    T = MaskedField{F}
#    A = typeof(field_array)
#    X = typeof(_axes)
#    new{T,N,A,X}(_axes,field_array,mask)
#  end
#end
#
#Base.size(a::MaskedFieldArray) = map(length,Base.axes(a))
#Base.axes(a::MaskedFieldArray) = a.axes
#Base.IndexStyle(::Type{<:MaskedFieldArray}) = IndexCartesian()
#function Base.getindex(a::MaskedFieldArray{T,N},i::Vararg{Integer,N}) where {T,N}
#  if a.mask
#    MaskedField(testitem(a.field_array),a.mask)
#  else
#    MaskedField(a.field_array[i...],a.mask)
#  end
#end
#
##function return_cache(a::MaskedFieldArray,x::Point)
##  cf = return_cache(a.field_array,x)
##  fx = return_value(a.field_array,x)
##  cz = CachedArray(similar(fx,eltype(fx),a.axes))
##  (cf,cz)
##end
##
##@inline function evaluate!(cache,a::MaskedFieldArray,x::Point)
##  cf, cz = cache
##  if a.mask
##    setaxes!(cz,a.axes)
##    r = cz.array
##    fill_entries!(r,zero(eltype(r)))
##    r
##  else
##    evaluate!(cf,a.field_array,x)
##  end
##end
##
##function return_cache(a::MaskedFieldArray,x::AbstractVector{<:Point})
##  cf = return_cache(a.field_array,x)
##  fx = return_value(a.field_array,x)
##  rx = similar_range(first(a.axes),length(x))
##  shape = (,)
##  cz = CachedArray(similar(fx,eltype(fx),a.axes))
##  (cf,cz)
##end
##
##@inline function evaluate!(cache,a::MaskedFieldArray,x::AbstractVector{<:Point})
##  cf, cz = cache
##  if a.mask
##    setaxes!(cz,a.axes)
##    r = cz.array
##    fill_entries!(r,zero(eltype(r)))
##    r
##  else
##    evaluate!(cf,a.field_array,x)
##  end
##end

#function _x_range(x:AbstractVector{<:Point},ran::Tuple{Base.OneTo})
#  Base.OneTo(length(x))
#end
#
#function _x_range(x:AbstractVector{<:Point},ran::Tuple{BlockedUnitRange})
#  blockedrange([length(x)])
#end
#
#function _x_range(x:AbstractVector{<:Point},ran::Tuple{BlockedUnitRange})
#_x_range(x:AbstractVector{<:Point},ran::Tuple{BlockedUnitRange})
#  blockedrange([length(x)])
#end
#  blockedrange([length(x)])
#end
#
#function _axes_with_x(x:AbstractVector{<:Point},ran::NTuple{N,<:BlockedUnitRange} where N)
#  np = length(x)
#  (blockedrange([np]),ran...)
#end
#
#function _new_axes(x,ran::NTuple{N,<:MultiLevelBlockedUnitRange} where N)
#  np = length(x)
#  a = _new_axes(x,first())
#  r = blockedrange([np])
#  (blockedrange([r]),ran...)
#end


#
#@inline function evaluate!(cache,a::BlockFieldArray,x::Point)
#  fc, blockids, cr = cache
#  fx = evaluate!(fc,a.field_array,x)
#  blockids[1] = a.blockid
#  evaluate!(cr,BlockVectorCoo,a.axes,blockids,fx)
#end


#struct BlockFieldArray{T,N,A,X} <: AbstractBlockArray{T,N}
#  axes::X
#  blockid::NTuple{N,Int}
#  field_array::A
#  function BlockFieldArray(
#    axes::NTuple{N},
#    blockid::NTuple{N,Int},
#    field_array::AbstractArray{T,N}) where {T<:Field,N}
#
#    @check begin 
#      msg = "The given field_array and axes are incompatible."
#      blocks_equal(axes(field_array),map(local_range,axes,blockid)) msg
#    end
#
#    A = typeof(field_array)
#    X = typeof(axes)
#    new{T,N,A,X}(axes,blockid,field_array)
#  end
#end
#
#Base.size(a::BlockFieldArray) = map(length,Base.axes(a))
#Base.axes(a::BlockFieldArray) = a.axes
#Base.IndexStyle(::Type{<:BlockFieldArray}) = IndexCartesian()
#Base.getindex(a::BlockFieldArray{T,N},i::Vararg{Integer,N}) where {T,N} = @notimplemented
#Base.setindex!(a::BlockFieldArray{T,N},v,i::Vararg{Integer,N}) where {T,N} = @notimplemented
#
#is_zero_block(a::BlockFieldArray{T,N},i::Vararg{Integer,N}) where {T,N} = i!=a.blockid
#BlockArrays.eachblock(a::BlockFieldArray) = ( a[I]  for I in eachblockid(a) )
#function BlockArrays.getblock(a::BlockFieldArray{T,N}, block::Vararg{Integer, N}) where {T,N}
#  if block == a.blockid
#    a.field_array
#  else
#    ai = testitem(a.field_array)
#    laxes = 
#    fill(zero(a),)
#
#  end
#end
#
#function return_cache(a::BlockFieldArray,x::Point)
#  fc = return_cache(a.field_array,x)
#  fx = return_value(a.field_array,x)
#  blockids = [a.blockid]
#  cr = return_cache(BlockVectorCoo,a.axes,blockids,fx)
#  (fc,blockids,cr)
#end
#
#@inline function evaluate!(cache,a::BlockFieldArray,x::Point)
#  fc, blockids, cr = cache
#  fx = evaluate!(fc,a.field_array,x)
#  blockids[1] = a.blockid
#  evaluate!(cr,BlockVectorCoo,a.axes,blockids,fx)
#end
