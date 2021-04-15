
struct MultiFieldFEBasis{BS,DS} <: FESpaces.FEBasis
  field::Int
  basis::SingleFieldFEBasis{BS,DS}
  spaces::Vector{<:SingleFieldFESpace}
end

Geometry.get_triangulation(f::MultiFieldFEBasis) = get_triangulation(f.basis)
FESpaces.BasisStyle(::Type{MultiFieldFEBasis{BS,DS}}) where {BS,DS} = BS()
CellData.DomainStyle(::Type{MultiFieldFEBasis{BS,DS}}) where {BS,DS} = DS()
Fields.gradient(a::MultiFieldFEBasis) =  MultiFieldFEBasis(a.field,∇(a.basis),a.spaces)
Fields.∇∇(a::MultiFieldFEBasis) =  MultiFieldFEBasis(a.field,∇∇(a.basis),a.spaces)

function _get_cell_axes(f::MultiFieldFEBasis,trian::Triangulation)
  nfields = length(f.spaces)
  fields =1:nfields
  field_i = f.field
  V_i = f.spaces[field_i]
  cell_axes = map(fields) do field_j
    V_j = f.spaces[field_j]
    trian_j = get_triangulation(V_j)
    if have_compatible_domains(trian_j,trian) ||
      have_compatible_domains(
      trian_j,get_background_triangulation(trian)) ||
      have_compatible_domains(
        trian_j,get_background_triangulation(get_background_triangulation(trian))) ||
      Geometry.is_included(trian,trian_j)
      cell_dofs_j = get_cell_dof_ids(V_j,trian)
      lazy_map(axes,cell_dofs_j)
    else
      cell_dofs_i = get_cell_dof_ids(V_i)
      dofs_i = testitem(cell_dofs_i)
      axes_i = axes(dofs_i)
      axes_i_empty = (similar_range(first(axes_i),0),)
      Fill(axes_i_empty,num_cells(trian))
    end
  end
  lazy_map(_multifield_axes_dofs,cell_axes...)
end

function CellData.get_data(f::MultiFieldFEBasis)
  nfields = length(f.spaces)
  field_i = f.field
  V_i = f.spaces[field_i]
  trian_i = get_triangulation(V_i)
  cell_axes_test = _get_cell_axes(f,trian_i)
  dv = f.basis
  if BasisStyle(f) == TestBasis()
    bsize = (nfields,)
    blockids = [(field_i,)]
    cell_axes = cell_axes_test
  elseif BasisStyle(f) == TrialBasis()
    bsize = (1,nfields)
    blockids = [(1,field_i)]
    cell_axes = lazy_map(cell_axes_test) do axs
      r = similar_range(first(axs),1)
      (r,axs...)
    end
  else
    @unreachable
  end
  cell_basis = lazy_map(Fields.BlockFieldArrayCooMap(bsize,blockids),cell_axes,get_data(dv))
end

function CellData.change_domain(
  f::MultiFieldFEBasis,trian::Triangulation,domain::DomainStyle)
  if have_compatible_domains(get_triangulation(f),trian)
    if DomainStyle(f) == domain
      return change_domain(f,domain)
    else
      @notimplemented
    end
  end
  nfields = length(f.spaces)
  field_i = f.field
  V_i = f.spaces[field_i]
  trian_i = get_triangulation(V_i)
  cell_axes_test = _get_cell_axes(f,trian)
  bsize = (nfields,)
  dv = CellData.change_domain(f.basis,trian,domain)
  if BasisStyle(f) == TestBasis()
    bsize = (nfields,)
    blockids = [(field_i,)]
    cell_axes = cell_axes_test
  elseif BasisStyle(f) == TrialBasis()
    bsize = (1,nfields)
    blockids = [(1,field_i)]
    cell_axes = lazy_map(cell_axes_test) do axs
      r = similar_range(first(axs),1)
      (r,axs...)
    end
  else
    @unreachable
  end
  cell_basis = lazy_map(Fields.BlockFieldArrayCooMap(bsize,blockids),cell_axes,get_data(dv))
  SingleFieldFEBasis(cell_basis,get_triangulation(dv),BasisStyle(f),DomainStyle(dv))
end

function CellData.change_domain(
  f::MultiFieldFEBasis,trian::SkeletonTriangulation,target_domain::DomainStyle)
  @unreachable """\n
  It is not possible to use the given MultiFieldFEBasis on a SkeletonTriangulation.
  Make sure that you are specifying which of the two possible traces,
  either plus (aka ⁺) or minus (aka ⁻) you want to use.
  """
end

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

function FESpaces.get_free_dof_ids(f::MultiFieldFESpace)
  block_dof_ids = AbstractUnitRange[]
  for U in f.spaces
    push!(block_dof_ids,get_free_dof_ids(U))
  end
  MultiLevelBlockedUnitRange(block_dof_ids)
end

FESpaces.get_dof_value_type(f::MultiFieldFESpace{MS,CS,V}) where {MS,CS,V} = eltype(V)

FESpaces.get_vector_type(f::MultiFieldFESpace) = f.vector_type

FESpaces.ConstraintStyle(::Type{MultiFieldFESpace{S,B,V}}) where {S,B,V} = B()

function FESpaces.get_cell_shapefuns(f::MultiFieldFESpace)
  field_to_febasis = map(get_cell_shapefuns,f.spaces)
  nfields = length(f.spaces)
  fields =1:nfields
  all_febases = map(fields) do field_i
    MultiFieldFEBasis(field_i,field_to_febasis[field_i],f.spaces)
  end
  MultiFieldCellField(all_febases)
end

function _multifield_axes_dofs(axs...)
  rs = map(first,axs)
  (append_ranges([rs...]),)
end

function FESpaces.get_cell_shapefuns_trial(f::MultiFieldFESpace)
  field_to_febasis = map(get_cell_shapefuns_trial,f.spaces)
  nfields = length(f.spaces)
  fields =1:nfields
  all_febases = map(fields) do field_i
    MultiFieldFEBasis(field_i,field_to_febasis[field_i],f.spaces)
  end
  MultiFieldCellField(all_febases)
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
    ma = cell_values.maps.value
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
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get
  the constrain flag for the corresponding cells.
  """
  trians = map(get_triangulation,f.spaces)
  trian = first(trians)
  @check all(map(t->have_compatible_domains(t,trian),trians)) msg
  get_cell_isconstrained(f,trian)
end

function FESpaces.get_cell_isconstrained(f::MultiFieldFESpace,trian::Triangulation)
  data = map(f.spaces) do space
    trian_i = get_triangulation(space)
    if have_compatible_domains(trian_i,trian) ||
      have_compatible_domains(trian_i,get_background_triangulation(trian)) ||
      have_compatible_domains(
        trian_i,get_background_triangulation(get_background_triangulation(trian))) ||
      Geometry.is_included(trian,trian_i)
      get_cell_isconstrained(space,trian)
    else
      Fill(false,num_cells(trian))
    end
  end
  lazy_map( (args...) -> +(args...)>0,  data...)
end

function FESpaces.get_cell_isconstrained(f::MultiFieldFESpace,trian::SkeletonTriangulation)
  plus = get_cell_isconstrained(f,trian.plus)
  minus = get_cell_isconstrained(f,trian.minus)
  lazy_map((l,r)-> l||r,plus,minus)
end

function FESpaces.get_cell_is_dirichlet(f::MultiFieldFESpace)
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get
  the constrain flag for the corresponding cells.
  """
  trians = map(get_triangulation,f.spaces)
  trian = first(trians)
  @check all(map(t->have_compatible_domains(t,trian),trians)) msg
  get_cell_is_dirichlet(f,trian)
end

function FESpaces.get_cell_is_dirichlet(f::MultiFieldFESpace,trian::Triangulation)
  data = map(f.spaces) do space
    trian_i = get_triangulation(space)
    if have_compatible_domains(trian_i,trian) ||
      have_compatible_domains(trian_i,get_background_triangulation(trian)) ||
      have_compatible_domains(
        trian_i,get_background_triangulation(get_background_triangulation(trian))) ||
      Geometry.is_included(trian,trian_i)
      get_cell_is_dirichlet(space,trian)
    else
      Fill(false,num_cells(trian))
    end
  end
  lazy_map( (args...) -> +(args...)>0,  data...)
end

function FESpaces.get_cell_is_dirichlet(f::MultiFieldFESpace,trian::SkeletonTriangulation)
  plus = get_cell_is_dirichlet(f,trian.plus)
  minus = get_cell_is_dirichlet(f,trian.minus)
  lazy_map((l,r)-> l||r,plus,minus)
end

function FESpaces.get_cell_constraints(f::MultiFieldFESpace)
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get
  the constrains for the corresponding cells.
  """
  trians = map(get_triangulation,f.spaces)
  trian = first(trians)
  @check all(map(t->have_compatible_domains(t,trian),trians)) msg
  get_cell_constraints(f,trian)
end

function FESpaces.get_cell_constraints(f::MultiFieldFESpace,trian::Triangulation)
  blocks = map(f.spaces) do space
    trian_i = get_triangulation(space)
    if have_compatible_domains(trian_i,trian) ||
      have_compatible_domains(trian_i,get_background_triangulation(trian)) ||
      have_compatible_domains(
        trian_i,get_background_triangulation(get_background_triangulation(trian))) ||
      Geometry.is_included(trian,trian_i)
      cell_constrs = get_cell_constraints(space,trian)
    else
      cell_constrs_i = get_cell_constraints(space)
      T = eltype(eltype(cell_constrs_i))
      cell_constrs = Fill(Matrix{T}(undef,0,0),num_cells(trian))
    end
    cell_constrs
  end
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

function FESpaces.get_cell_constraints(f::MultiFieldFESpace,trian::SkeletonTriangulation)
  cell_constraints_plus = get_cell_constraints(f,trian.plus)
  cell_constraints_minus = get_cell_constraints(f,trian.minus)
  cell_axes_plus = lazy_map(axes,cell_constraints_plus)
  cell_axes_minus = lazy_map(axes,cell_constraints_minus)
  cell_axes = lazy_map(cell_axes_plus,cell_axes_minus) do axp, axm
    r1 = append_ranges([axp[1],axm[1]])
    r2 = append_ranges([axp[2],axm[2]])
    (r1,r2)
  end
  lazy_map(
    BlockArrayCooMap((2,2),[(1,1),(2,2)]),
    cell_axes,
    cell_constraints_plus,
    cell_constraints_minus)
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace)
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get the DOF ids
  on top of the corresponding cells.
  """
  trians = map(get_triangulation,f.spaces)
  trian = first(trians)
  @check all(map(t->have_compatible_domains(t,trian),trians)) msg
  get_cell_dof_ids(f,trian)
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,trian::Triangulation)
  get_cell_dof_ids(f,trian,MultiFieldStyle(f))
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,trian::SkeletonTriangulation)
  cell_ids_plus = get_cell_dof_ids(f,trian.plus)
  cell_ids_minus = get_cell_dof_ids(f,trian.minus)
  cell_axes_plus = lazy_map(axes,cell_ids_plus)
  cell_axes_minus = lazy_map(axes,cell_ids_minus)
  cell_axes = lazy_map(cell_axes_plus,cell_axes_minus) do axp, axm
    (append_ranges([axp[1],axm[1]]),)
  end
  lazy_map(BlockArrayCooMap((2,),[(1,),(2,)]),cell_axes,cell_ids_plus,cell_ids_minus)
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,::Triangulation,::MultiFieldStyle)
  @notimplemented
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,trian::Triangulation,::ConsecutiveMultiFieldStyle)
  offsets = compute_field_offsets(f)
  spaces = f.spaces
  function fun(i,space)
    trian_i = get_triangulation(space)
    if have_compatible_domains(trian_i,trian) ||
      have_compatible_domains(trian_i,get_background_triangulation(trian)) ||
      have_compatible_domains(
        trian_i,get_background_triangulation(get_background_triangulation(trian))) ||
      Geometry.is_included(trian,trian_i)
      cell_dofs = get_cell_dof_ids(space,trian)
    else
      cell_dofs_i = get_cell_dof_ids(space)
      T = eltype(eltype(cell_dofs_i))
      cell_dofs = Fill(T[],num_cells(trian))
    end
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
