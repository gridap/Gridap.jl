
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
struct MultiFieldFESpace{MS<:MultiFieldStyle,CS<:ConstraintStyle,V} <: FESpace
  vector_type::Type{V}
  spaces::Vector{<:SingleFieldFESpace}
  multi_field_style::MS
  constraint_style::CS
  function MultiFieldFESpace(
    ::Type{V},
    spaces::Vector{<:SingleFieldFESpace},
    multi_field_style::MultiFieldStyle) where V
    @assert length(spaces) > 0

    MS = typeof(multi_field_style)
    if any( map(has_constraints,spaces) )
      constraint_style = Constrained()
    else
      constraint_style = UnConstrained()
    end
    CS = typeof(constraint_style)
    new{MS,CS,V}(V,spaces,multi_field_style,constraint_style)
  end
end

"""
    MultiFieldFESpace(spaces::Vector{<:SingleFieldFESpace})
"""
function MultiFieldFESpace(spaces::Vector{<:SingleFieldFESpace})
  Ts = map(get_dof_value_type,spaces)
  T = typeof(*(map(zero,Ts)...))
  MultiFieldFESpace(Vector{T},spaces,ConsecutiveMultiFieldStyle())
end

function MultiFieldFESpace(::Type{V},spaces::Vector{<:SingleFieldFESpace}) where V
  MultiFieldFESpace(V,spaces,ConsecutiveMultiFieldStyle())
end

MultiFieldStyle(::Type{MultiFieldFESpace{S,B,V}}) where {S,B,V} = S()
MultiFieldStyle(f::MultiFieldFESpace) = MultiFieldStyle(typeof(f))

# Implementation of FESpace

function FESpaces.get_triangulation(f::MultiFieldFESpace)
  s1 = first(f.spaces)
  trian = get_triangulation(s1)
  @check all(map(i->trian===get_triangulation(i),f.spaces))
  trian
end

function FESpaces.num_free_dofs(f::MultiFieldFESpace)
  n = 0
  for U in f.spaces
    n += num_free_dofs(U)
  end
  n
end

FESpaces.get_dof_value_type(f::MultiFieldFESpace{MS,CS,V}) where {MS,CS,V} = eltype(V)

FESpaces.get_vector_type(f::MultiFieldFESpace) = f.vector_type

FESpaces.ConstraintStyle(::Type{MultiFieldFESpace{S,B,V}}) where {S,B,V} = B()

function FESpaces.get_cell_shapefuns(f::MultiFieldFESpace)
 
  # Compute the cell axes
  blocks = map(get_cell_shapefuns,f.spaces)
  blocks_data = map(get_cell_data,blocks)
  cell_axes_blocks = map(i->lazy_map(axes,i),blocks_data)
  cell_axes = lazy_map(_multifield_axes_dofs,cell_axes_blocks...)

  # Compute the underlying single-fields
  nfields = length(f.spaces)
  all_bases = map(1:nfields) do i
    bsize = (nfields,)
    dv = blocks[i]
    cell_basis = lazy_map(Fields.BlockFieldArrayCooMap(bsize,[(i,)]),cell_axes,get_cell_data(dv))
    FEBasis(cell_basis,get_triangulation(dv),TestBasis(),DomainStyle(dv))
  end
  MultiFieldCellField(all_bases)
end

function _multifield_axes_dofs(axs...)
  rs = map(first,axs)
  (append_ranges([rs...]),)
end

function FESpaces.get_cell_shapefuns_trial(f::MultiFieldFESpace)
  blocks = get_cell_shapefuns(f).single_fields
  nfields = length(blocks)
  all_bases = map(1:nfields) do i
    dv = blocks[i]
    cell_basis = lazy_map(transpose,get_cell_data(dv))
    FEBasis(cell_basis,get_triangulation(dv),TrialBasis(),DomainStyle(dv))
  end
  MultiFieldCellField(all_bases)
end

function FESpaces.FEFunction(fe::MultiFieldFESpace, free_values)
  blocks = map(1:length(fe.spaces)) do i
    free_values_i = restrict_to_field(fe,free_values,i)
    FEFunction(fe.spaces[i],free_values_i)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

function FESpaces.EvaluationFunction(fe::MultiFieldFESpace, free_values)
  blocks = map(1:length(fe.spaces)) do i
    free_values_i = restrict_to_field(fe,free_values,i)
    EvaluationFunction(fe.spaces[i],free_values_i)
  end
  MultiFieldFEFunction(free_values,fe,blocks)
end

function CellData.CellField(fe::MultiFieldFESpace,cell_values)
  single_fields = map(1:length(fe.spaces)) do i
    cell_values_field = lazy_map(a->a[Block(i)],cell_values)
    CellField(fe.spaces[i],cell_values_field)
  end
  MultiFieldCellField(single_fields)
end

function CellData.CellField(
  fe::MultiFieldFESpace,cell_values::LazyArray{<:Fill{<:BlockArrayCooMap}})
  single_fields = map(1:length(fe.spaces)) do i
    ma = cell_values.g.value
    cell_values_field = Fields._get_cell_block(ma,cell_values,Block(i))
    CellField(fe.spaces[i],cell_values_field)
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

function FESpaces.get_cell_isconstrained(f::MultiFieldFESpace)
  data = map(get_cell_isconstrained,f.spaces)
  lazy_map( (args...) -> +(args...)>0,  data...)
end

function FESpaces.get_cell_constraints(f::MultiFieldFESpace)
  blocks = map(get_cell_constraints,f.spaces)
  blockids = [ (i,i) for i in 1:length(f.spaces)]
  nfields = length(f.spaces)
  bsize = (nfields,nfields)

  blocks_axes = map( i->lazy_map(axes,i), blocks)

  function multifield_axes_constraints(axs...)
    rs1 = map(j->j[1],axs)
    rs2 = map(j->j[2],axs)
    r1 = append_ranges([rs1...])
    r2 = append_ranges([rs2...])
    (r1,r2)
  end
  cell_axes = lazy_map(multifield_axes_constraints,blocks_axes...)
  lazy_map(BlockArrayCooMap(bsize,blockids),cell_axes,blocks...)
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace)
  get_cell_dof_ids(f,MultiFieldStyle(f))
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,::MultiFieldStyle)
  @notimplemented
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,::ConsecutiveMultiFieldStyle)
  offsets = compute_field_offsets(f)
  spaces = f.spaces
  function fun(i,space)
    cell_dofs = get_cell_dof_ids(space)
    if i == 1
      return cell_dofs
    end
    offset = Int32(offsets[i])
    o = Fill(offset,length(cell_dofs))
    lazy_map(Broadcasting(_sum_if_first_positive),cell_dofs,o)
  end
  blocks = [ fun(i,space) for (i,space) in enumerate(spaces) ]
  bids = [ (i,) for i in 1:length(spaces)]
  bsize = (length(f.spaces),)
  cell_axes_blocks = map(i->lazy_map(axes,i),blocks)
  cell_axes = lazy_map(_multifield_axes_dofs,cell_axes_blocks...)
  lazy_map(BlockArrayCooMap(bsize,bids),cell_axes,blocks...)
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
