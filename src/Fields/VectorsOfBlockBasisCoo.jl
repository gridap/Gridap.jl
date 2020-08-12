
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
  axs
end

Base.size(v::VectorOfBlockBasisCoo) = (length(first(v.blocks)) ,)
Base.getindex(v::VectorOfBlockBasisCoo,i::Integer) = BlockBasisCoo()

function evaluate_field_array(v::VectorOfBlockBasisCoo,x::AbstractArray)
  blocks = evaluate_field_arrays(v.blocks,x)
  VectorOfBlockArrayCoo(blocks,v.blockids,v.axs)
end

function field_array_gradient(v::VectorOfBlockBasisCoo)
  blocks = map(field_array_gradient,v.blocks)
  VectorOfBlockBasisCoo(blocks,v.blockids,v.axs)
end

