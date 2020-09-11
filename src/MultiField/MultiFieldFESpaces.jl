
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
    B = any( map(has_constraints,spaces) )
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

FESpaces.constraint_style(::Type{MultiFieldFESpace{S,B}}) where {S,B} = Val{B}()

function FESpaces.num_free_dofs(f::MultiFieldFESpace)
  n = 0
  for U in f.spaces
    n += num_free_dofs(U)
  end
  n
end

function FESpaces.get_cell_axes(f::MultiFieldFESpace)
  all_axes = map(get_cell_axes,f.spaces)
  create_array_of_blocked_axes(all_axes...)
end

function FESpaces.get_cell_axes_with_constraints(f::MultiFieldFESpace)
  all_axes = map(get_cell_axes_with_constraints,f.spaces)
  create_array_of_blocked_axes(all_axes...)
end

function FESpaces.get_cell_basis(f::MultiFieldFESpace)
  all_axes = map(i->get_cell_axes(get_cell_basis(i)),f.spaces)
  all_bases = [
    _insert_cell_basis_in_block(i,get_cell_basis(f.spaces[i]),all_axes)
    for i in 1:length(f.spaces) ]
  MultiFieldCellField(all_bases)
end

function _insert_cell_basis_in_block(i::Int,_cb,all_axes)
  cb = get_wrapped_cell_field(_cb)
  array =  insert_array_of_bases_in_block(i,get_array(cb),all_axes...)
  c = GenericCellField(array,array.axes,MetaSizeStyle(cb))
  wrap_cell_field(_cb,c)
end

function FESpaces.FEFunction(fe::MultiFieldFESpace, free_values)
  blocks = SingleFieldFEFunction[]
  for (field, U) in enumerate(fe.spaces)
    free_values_i = restrict_to_field(fe,free_values,field)
    uhi = FEFunction(U,free_values_i)
    push!(blocks,uhi)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

function FESpaces.EvaluationFunction(fe::MultiFieldFESpace, free_values)
  blocks = SingleFieldFEFunction[]
  for (field, U) in enumerate(fe.spaces)
    free_values_i = restrict_to_field(fe,free_values,field)
    uhi = EvaluationFunction(U,free_values_i)
    push!(blocks,uhi)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

function CellData.CellField(fe::MultiFieldFESpace,cell_values)
  single_fields = CellField[]
  for field in 1:num_fields(fe)
    cell_values_field = apply(a->a[Block(field)],cell_values)
    cf = CellField(fe.spaces[field],cell_values_field)
    push!(single_fields,cf)
  end
  MultiFieldCellField(single_fields)
end

function CellData.CellField(fe::MultiFieldFESpace,cell_values::VectorOfBlockArrayCoo)
  single_fields = CellField[]
  for field in 1:num_fields(fe)
    cf = CellField(fe.spaces[field],cell_values[Block(field)])
    push!(single_fields,cf)
  end
  MultiFieldCellField(single_fields)
end

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

function FESpaces.zero_free_values(fs::MultiFieldFESpace)
  zeros(num_free_dofs(fs))
end

function FESpaces.get_cell_isconstrained(f::MultiFieldFESpace)
  data = map(get_cell_isconstrained,f.spaces)
  apply( (args...) -> +(args...)>0,  data...)
end

function FESpaces.get_cell_constraints(f::MultiFieldFESpace)
  blocks = Tuple(map(get_cell_constraints,f.spaces))
  blockids = [ (i,i) for i in 1:length(f.spaces)]
  all_axs_rows = map(get_cell_axes_with_constraints,f.spaces)
  all_axs_cols = map(get_cell_axes,f.spaces)
  axs_rows = create_array_of_blocked_axes(all_axs_rows...)
  axs_cols = create_array_of_blocked_axes(all_axs_cols...)
  axs = apply((r,c) -> (r[1],c[1]),axs_rows,axs_cols)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

function FESpaces.get_cell_dofs(f::MultiFieldFESpace)
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
  blocks = Tuple([ fun(i,space) for (i,space) in enumerate(spaces) ])
  blockids = [ (i,) for i in 1:length(spaces)]
  axs = get_cell_axes(f)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

function _sum_if_first_positive(a,b)
  if a>0
    return a+b
  else
    return a
  end
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


# API for the ConsecutiveMultiFieldStyle
import Gridap.FESpaces: interpolate
import Gridap.FESpaces: interpolate_everywhere
import Gridap.FESpaces: interpolate_dirichlet

@deprecate(
  interpolate(fs::MultiFieldFESpace, object),
  interpolate(object, fs::MultiFieldFESpace)
)

@deprecate(
  interpolate_everywhere(fs::MultiFieldFESpace, object),
  interpolate_everywhere(object, fs::MultiFieldFESpace)
)

@deprecate(
  interpolate_dirichlet(fs::MultiFieldFESpace, object),
  interpolate_dirichlet(object, fs::MultiFieldFESpace)
)

"""
The resulting MultiFieldFEFunction is in the space (in particular it fulfills Dirichlet BCs
even in the case that the given cell field does not fulfill them)
"""
function FESpaces.interpolate(objects, fe::MultiFieldFESpace)
  free_values = zero_free_values(fe)
  blocks = SingleFieldFEFunction[]
  for (field, (U,object)) in enumerate(zip(fe.spaces,objects))
    free_values_i = restrict_to_field(fe,free_values,field)
    uhi = interpolate!(object, free_values_i,U)
    push!(blocks,uhi)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

"""
like interpolate, but also compute new degrees of freedom for the dirichlet component.
The resulting MultiFieldFEFunction does not necessary belongs to the underlying space
"""
function FESpaces.interpolate_everywhere(objects, fe::MultiFieldFESpace)
  free_values = zero_free_values(fe)
  blocks = SingleFieldFEFunction[]
  for (field, (U,object)) in enumerate(zip(fe.spaces,objects))
    free_values_i = restrict_to_field(fe,free_values,field)
    dirichlet_values_i = zero_dirichlet_values(U)
    uhi = interpolate_everywhere!(object, free_values_i,dirichlet_values_i,U)
    push!(blocks,uhi)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

"""
"""
function FESpaces.interpolate_dirichlet(objects, fe::MultiFieldFESpace)
  free_values = zero_free_values(fe)
  blocks = SingleFieldFEFunction[]
  for (field, (U,object)) in enumerate(zip(fe.spaces,objects))
    free_values_i = restrict_to_field(fe,free_values,field)
    dirichlet_values_i = zero_dirichlet_values(U)
    uhi = interpolate_dirichlet!(object, free_values_i,dirichlet_values_i,U)
    push!(blocks,uhi)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

