
struct BlockBasisCoo <: Field end

function evaluate_field!(cache,b::BlockBasisCoo,x)
  msg =
  """
  Not implemented! We not support iteration of vectors of block bases before
  evaluation. This is unlikely to be be needed in the future.
  """
  @notimplemented msg
end

struct VectorOfBlockBasisCoo <: AbstractVector{BlockBasisCoo}
  blocks::Tuple
  blockids
  axes
end

Base.size(v::VectorOfBlockBasisCoo) = (length(first(v.blocks)) ,)
Base.getindex(v::VectorOfBlockBasisCoo,i::Integer) = BlockBasisCoo()

function evaluate_field_array(v::VectorOfBlockBasisCoo,x::AbstractArray)
  blocks = evaluate_field_arrays(v.blocks,x)
  axs = apply(_new_axes,x,v.axes)
  blockids = map(i->(1,i...),v.blockids)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

function field_array_gradient(v::VectorOfBlockBasisCoo)
  blocks = map(field_array_gradient,v.blocks)
  VectorOfBlockBasisCoo(blocks,v.blockids,v.axes)
end

function reindex(v::VectorOfBlockBasisCoo,j_to_i::AbstractArray)
  blocks = map(b->reindex(b,j_to_i),v.blocks)
  axs = reindex(v.axes,j_to_i)
  VectorOfBlockBasisCoo(blocks,v.blockids,axs)
end

function reindex(v::VectorOfBlockBasisCoo,j_to_i::IdentityVector)
  @assert length(v) == length(j_to_i)
  v
end

function compose_field_arrays(v::VectorOfBlockBasisCoo,f)
  blocks = map(b->compose_field_arrays(b,f),v.blocks)
  VectorOfBlockBasisCoo(blocks,v.blockids,v.axes)
end

function insert_array_of_bases_in_block(i::Integer,a,ax1,ax2...)
  blocks = (a,)
  blockids = _compute_blockids(eltype(ax1),i)
  axs = create_array_of_blocked_axes(ax1,ax2...)
  VectorOfBlockBasisCoo(blocks,blockids,axs)
end

#TODO perhaps this should go into Arrays/BlockArraysCoo.jl
function create_array_of_blocked_axes(axs...)
  axs = apply(_cat_axes,axs...)
end

function _compute_blockids(::Type{<:NTuple{1}},i)
  [(i,)]
end

function _compute_blockids(::Type{<:NTuple{2}},i)
  [(1,i)]
end

function _cat_axes(a::NTuple{1},b::NTuple{1})
  (_blockedrange([a[1],b[1]]),)
end

function _cat_axes(a::NTuple{1}...)
  (_blockedrange([map(i->i[1],a)...]),)
end

function _cat_axes(a::NTuple{2},b::NTuple{2})
  @assert length(a[1]) == 1
  @assert length(b[1]) == 1
  ran = (_blockedrange([a[2],b[2]]),)
  _add_singleton_block(ran)
end

function _cat_axes(a::NTuple{2}...)
  ran = (_blockedrange([map(i->i[2],a)...]),)
  _add_singleton_block(ran)
end

function _new_axes(x,ran::NTuple{N,<:BlockedUnitRange} where N)
  np = length(x)
  (blockedrange([np]),ran...)
end

function _new_axes(x,ran::NTuple{N,<:TwoLevelBlockedUnitRange} where N)
  np = length(x)
  r = blockedrange([np])
  (blockedrange([r]),ran...)
end

function _add_singleton_block(ran::NTuple{N,<:Base.OneTo} where N)
  (Base.OneTo(1),ran...)
end

function _add_singleton_block(ran::NTuple{N,<:BlockedUnitRange} where N)
  (blockedrange([1]),ran...)
end

function _add_singleton_block(ran::NTuple{N,<:TwoLevelBlockedUnitRange} where N)
  r = blockedrange([1])
  (blockedrange([r]),ran...)
end

function _blockedrange(a::Vector{<:Base.OneTo})
  blockedrange(map(length,a))
end

function _blockedrange(a::Vector{<:BlockedUnitRange})
  blockedrange(a)
end


