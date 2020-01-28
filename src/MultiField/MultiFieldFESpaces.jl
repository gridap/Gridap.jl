
abstract type MultiFieldStyle end

struct ConsequtiveMultiFieldStyle <: MultiFieldStyle end

struct StridedMultiFieldStyle <: MultiFieldStyle end

struct MultiFieldFESpace{S<:MultiFieldStyle} <: FESpace
  spaces::Vector{<:SingleFieldFESpace}
  multi_field_style::S
end

MultiFieldStyle(::Type{MultiFieldFESpace{S}}) where S = S()

MultiFieldStyle(f::MultiFieldFESpace) = MultiFieldStyle(typeof(f))


# Implementation of FESpace

function num_free_dofs(f::MultiFieldFESpace)
  n = 0
  for U in f.spaces
    n += num_free_dofs(U)
  end
  n
end

function get_cell_basis(f::MultiFieldFESpace)
  blocks = [ get_cell_basis(U) for U in f.spaces ]
  MultiCellBasis(blocks)
end

function FEFunction(fe::MultiFieldFESpace, free_values)
  blocks = SingleFieldFEFunction[]
  for (field, U) in enumerate(fe.spaces)
    free_values_i = restrict_to_field(fe,free_values,field)
    uhi = FEFunction(U,free_values_i)
    push!(blocks,uhi)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

function zero_free_values(::Type{T},fs::MultiFieldFESpace) where T
  zeros(T,num_free_dofs(fs))
end

function apply_constraints_matrix_cols(f::MultiFieldFESpace,cellmat::MultiCellArray,cellids::AbstractVector)
  blocks = cellmat.blocks
  block_ids = cellmat.block_ids
  spaces = f.spaces
  function fun(i,block)
    field_id_rows, field_id_cols = block_ids[i]
    apply_constraints_matrix_cols(spaces[field_id_cols],block,cellids)
  end
  new_blocks = (  fun(i,block) for (i,block) in enumerate(blocks) )
  MultiCellArray(tuple(new_blocks...),block_ids)
end

function apply_constraints_matrix_rows(f::MultiFieldFESpace,cellmat::MultiCellArray,cellids::AbstractVector)
  blocks = cellmat.blocks
  block_ids = cellmat.block_ids
  spaces = f.spaces
  function fun(i,block)
    field_id_rows, field_id_cols = block_ids[i]
    apply_constraints_matrix_rows(spaces[field_id_rows],block,cellids)
  end
  new_blocks = (  fun(i,block) for (i,block) in enumerate(blocks) )
  MultiCellArray(tuple(new_blocks...),block_ids)
end

function apply_constraints_vector(f::MultiFieldFESpace,cellvec::MultiCellArray,cellids::AbstractVector)
  blocks = cellvec.blocks
  block_ids = cellvec.block_ids
  spaces = f.spaces
  function fun(i,block)
   field_id = block_ids[i]
   apply_constraints_vector(spaces[field_id],block,cellids)
  end
  new_blocks = (  fun(i,block) for (i,block) in enumerate(blocks) )
  MultiCellArray(tuple(new_blocks...),block_ids)
end

# API for multi field case

function num_fields(f::MultiFieldFESpace)
  length(f.spaces)
end

Base.iterate(m::MultiFieldFESpace) = iterate(m.spaces)

Base.iterate(m::MultiFieldFESpace,state) = iterate(m.spaces,state)

Base.getindex(m::MultiFieldFESpace,field_id::Integer) = m.spaces[field_id]

function restrict_to_field(f::MultiFieldFESpace,free_values::AbstractVector,field::Integer)
  _restrict_to_field(f,MultiFieldStyle(f),free_values,field)
end

function  _restrict_to_field(f,mf_style,free_values,field)
  @notimplemented
end

function  _restrict_to_field(f,::ConsequtiveMultiFieldStyle,free_values,field)
  offsets = compute_field_offsets(f)
  U = f.spaces
  pini = offsets[field] + 1
  pend = offsets[field] + num_free_dofs(U[field])
  SubVector(free_values,pini,pend)
end

# API for the ConsequtiveMultiFieldStyle

function compute_field_offsets(f::MultiFieldFESpace)
  @assert MultiFieldStyle(f) == ConsequtiveMultiFieldStyle()
  U = f.spaces
  n = length(U)
  offsets = zeros(Int,n)
  for i in 1:(n-1)
    Ui = U[i]
    offsets[i+1] = offsets[i] + num_free_dofs(Ui)
  end
  offsets
end

