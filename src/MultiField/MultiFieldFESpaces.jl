
abstract type MultiFieldStyle end

struct ConsecutiveMultiFieldStyle <: MultiFieldStyle end

struct StridedMultiFieldStyle <: MultiFieldStyle end

"""
    struct MultiFieldFESpace{S<:MultiFieldStyle,B} <: FESpace
      spaces::Vector{<:SingleFieldFESpace}
      multi_field_style::S
      constraint_style::Val{B}
    end
"""
struct MultiFieldFESpace{S<:MultiFieldStyle,B} <: FESpace
  spaces::Vector{<:SingleFieldFESpace}
  multi_field_style::S
  constraint_style::Val{B}
  function MultiFieldFESpace(spaces::Vector{<:SingleFieldFESpace},mfs::MultiFieldStyle)
    S = typeof(mfs)
    msg = "MultiFieldFESpace only implemented when the underlying fields have no constraints."
    @notimplementedif any( map(has_constraints,spaces) ) msg
    B = false
    cs = Val{B}()
    new{S,B}(spaces,mfs,cs)
  end
end

"""
    MultiFieldFESpace(spaces::Vector{<:SingleFieldFESpace})
"""
function MultiFieldFESpace(spaces::Vector{<:SingleFieldFESpace})
  MultiFieldFESpace(spaces,ConsecutiveMultiFieldStyle())
end

MultiFieldStyle(::Type{MultiFieldFESpace{S,B}}) where {S,B} = S()

MultiFieldStyle(f::MultiFieldFESpace) = MultiFieldStyle(typeof(f))

# Implementation of FESpace

constraint_style(::Type{MultiFieldFESpace{S,B}}) where {S,B} = Val{B}()

function num_free_dofs(f::MultiFieldFESpace)
  n = 0
  for U in f.spaces
    n += num_free_dofs(U)
  end
  n
end

function num_free_dofs(spaces::Vector{<:SingleFieldFESpace})
  f = MultiFieldFESpace(spaces)
  num_free_dofs(f)
end

function get_cell_basis(f::MultiFieldFESpace)
  blocks = [ get_cell_basis(U) for U in f.spaces ]
  MultiCellBasis(blocks)
end

function get_cell_basis(spaces::Vector{<:SingleFieldFESpace})
  f = MultiFieldFESpace(spaces)
  get_cell_basis(f)
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

function FEFunction(spaces::Vector{<:SingleFieldFESpace}, free_values)
  f = MultiFieldFESpace(spaces)
  FEFunction(f,free_values)
end

function EvaluationFunction(fe::MultiFieldFESpace, free_values)
  blocks = SingleFieldFEFunction[]
  for (field, U) in enumerate(fe.spaces)
    free_values_i = restrict_to_field(fe,free_values,field)
    uhi = EvaluationFunction(U,free_values_i)
    push!(blocks,uhi)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

function EvaluationFunction(spaces::Vector{<:SingleFieldFESpace}, free_values)
  f = MultiFieldFESpace(spaces)
  EvaluationFunction(f,free_values)
end

function zero_free_values(::Type{T},fs::MultiFieldFESpace) where T
  zeros(T,num_free_dofs(fs))
end

function zero_free_values(::Type{T},spaces::Vector{<:SingleFieldFESpace}) where T
  f = MultiFieldFESpace(spaces)
  zero_free_values(T,f)
end

# API for multi field case

"""
    num_fields(f::MultiFieldFESpace)
"""
function num_fields(f::MultiFieldFESpace)
  length(f.spaces)
end

Base.iterate(m::MultiFieldFESpace) = iterate(m.spaces)

Base.iterate(m::MultiFieldFESpace,state) = iterate(m.spaces,state)

Base.getindex(m::MultiFieldFESpace,field_id::Integer) = m.spaces[field_id]

"""
    restrict_to_field(f::MultiFieldFESpace,free_values::AbstractVector,field::Integer)
"""
function restrict_to_field(f::MultiFieldFESpace,free_values::AbstractVector,field::Integer)
  _restrict_to_field(f,MultiFieldStyle(f),free_values,field)
end

function  _restrict_to_field(f,::MultiFieldStyle,free_values,field)
  @notimplemented
end

function  _restrict_to_field(f,::ConsecutiveMultiFieldStyle,free_values,field)
  offsets = compute_field_offsets(f)
  U = f.spaces
  pini = offsets[field] + 1
  pend = offsets[field] + num_free_dofs(U[field])
  SubVector(free_values,pini,pend)
end

function get_cell_dofs(f::MultiFieldFESpace)
  _get_cell_dofs(f,MultiFieldStyle(f))
end

function _get_cell_dofs(f,::MultiFieldStyle)
  @notimplemented
end


function _get_cell_dofs(f,::ConsecutiveMultiFieldStyle)
  offsets = compute_field_offsets(f)
  spaces = f.spaces
  function fun(i,space)
    cell_dofs = get_cell_dofs(space)
    if i == 1
      return cell_dofs
    end
    offset = offsets[i]
    o = Fill(offset,length(cell_dofs))
    apply(elem(_sum_if_first_positive),cell_dofs,o)
  end
  blocks = [ fun(i,space) for (i,space) in enumerate(spaces) ]
  block_ids = [ (i,) for i in 1:length(spaces)]
  MultiFieldCellArray(tuple(blocks...),block_ids)
end

function _sum_if_first_positive(a,b)
  if a>0
    return a+b
  else
    return a
  end
end

# API for the ConsecutiveMultiFieldStyle

"""
    compute_field_offsets(f::MultiFieldFESpace)
"""
function compute_field_offsets(f::MultiFieldFESpace)
  @assert MultiFieldStyle(f) == ConsecutiveMultiFieldStyle()
  U = f.spaces
  n = length(U)
  offsets = zeros(Int,n)
  for i in 1:(n-1)
    Ui = U[i]
    offsets[i+1] = offsets[i] + num_free_dofs(Ui)
  end
  offsets
end

