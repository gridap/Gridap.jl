
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
  ranges
end

Base.size(v::VectorOfBlockBasisCoo) = (length(first(v.blocks)) ,)
Base.getindex(v::VectorOfBlockBasisCoo,i::Integer) = BlockBasisCoo()

function evaluate_field_array(v::VectorOfBlockBasisCoo,x::AbstractArray)
  blocks = evaluate_field_arrays(v.blocks,x)
  axs = apply(_new_axes,x,v.ranges)
  VectorOfBlockArrayCoo(blocks,v.blockids,axs)
end

function _new_axes(x,ran::BlockedUnitRange)
  np = length(x)
  (blockedrange([np]),ran)
end

function _new_axes(x,ran::TwoLevelBlockedUnitRange)
  np = length(x)
  r = blockedrange([np])
  (blockedrange([r]),ran)
end

function field_array_gradient(v::VectorOfBlockBasisCoo)
  blocks = map(field_array_gradient,v.blocks)
  VectorOfBlockBasisCoo(blocks,v.blockids,v.ranges)
end

