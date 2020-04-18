
"""
    struct MultiFieldFEFunction
      # private fields
    end
"""
struct MultiFieldFEFunction
  free_values::AbstractVector
  space::MultiFieldFESpace
  blocks::Vector{<:SingleFieldFEFunction}
  cell_vals::MultiFieldCellArray

  function MultiFieldFEFunction(
    free_values::AbstractVector,
    space::MultiFieldFESpace,
    blocks::Vector{<:SingleFieldFEFunction})

    _blocks = tuple(map(get_cell_values,blocks)...)
    _block_ids = [ (i,) for i in 1:length(blocks) ]
    cell_vals = MultiFieldCellArray(_blocks,_block_ids)

    new(free_values,space,blocks,cell_vals)
  end
end

function MultiFieldFEFunction(
  space::MultiFieldFESpace,
  blocks::Vector{<:SingleFieldFEFunction})
  fv = zero_free_values(space)
  xh0 = MultiFieldFEFunction(fv,X0,blocks)
end

FEFunctionStyle(::Type{MultiFieldFEFunction}) = Val{true}()

get_free_values(f::MultiFieldFEFunction) = f.free_values

get_fe_space(f::MultiFieldFEFunction) = f.space

"""
    num_fields(m::MultiFieldFEFunction)
"""
num_fields(m::MultiFieldFEFunction) = length(m.blocks)

Base.iterate(m::MultiFieldFEFunction) = iterate(m.blocks)

Base.iterate(m::MultiFieldFEFunction,state) = iterate(m.blocks,state)

Base.getindex(m::MultiFieldFEFunction,field_id::Integer) = m.blocks[field_id]

function restrict(a::MultiFieldFEFunction,trian::Triangulation)
  f = (ai) -> restrict(ai,trian)
  blocks = map(f,a.blocks)
  blocks
end

function get_cell_values(f::MultiFieldFEFunction)
  f.cell_vals
end
