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
  get_free_dof_ids(f,MultiFieldStyle(f))
end

function FESpaces.get_free_dof_ids(f::MultiFieldFESpace,::MultiFieldStyle)
  @abstractmethod
end

function FESpaces.get_free_dof_ids(f::MultiFieldFESpace,::ConsecutiveMultiFieldStyle)
  block_num_dofs = Int[]
  for U in f.spaces
    push!(block_num_dofs,num_free_dofs(U))
  end
  BlockArrays.blockedrange(block_num_dofs)
end

FESpaces.get_dof_value_type(f::MultiFieldFESpace{MS,CS,V}) where {MS,CS,V} = eltype(V)

FESpaces.get_vector_type(f::MultiFieldFESpace) = f.vector_type

FESpaces.ConstraintStyle(::Type{MultiFieldFESpace{S,B,V}}) where {S,B,V} = B()

struct MultiFieldFEBasisComponent{B} <: FEBasis
  cell_basis::AbstractArray
  single_field::B
  fieldid::Int
  nfields::Int
  function MultiFieldFEBasisComponent(
    single_field::SingleFieldFEBasis,fieldid::Integer,nfields::Integer)
    function block_dofs(cell_bs,::TestBasis,fieldid,nfields)
      cell_basis = lazy_map(BlockMap(nfields,fieldid),cell_bs)
    end
    function block_dofs(cell_bs,::TrialBasis,fieldid,nfields)
      cell_basis = lazy_map(BlockMap((1,nfields),fieldid),cell_bs)
    end
    B = typeof(single_field)
    cell_bs = get_data(single_field)
    bsty = BasisStyle(single_field)
    cell_basis = block_dofs(cell_bs,bsty,fieldid,nfields)
    new{B}(cell_basis,single_field,fieldid,nfields)
  end
end

CellData.get_data(f::MultiFieldFEBasisComponent) = f.cell_basis
CellData.get_triangulation(f::MultiFieldFEBasisComponent) = get_triangulation(f.single_field)
FESpaces.BasisStyle(::Type{<:MultiFieldFEBasisComponent{B}}) where B = BasisStyle(B)
CellData.DomainStyle(::Type{<:MultiFieldFEBasisComponent{B}}) where B = DomainStyle(B)
# Do not copy discretisation data
Base.copy(f::MultiFieldFEBasisComponent) = f
function FESpaces.CellData.similar_cell_field(
  f::MultiFieldFEBasisComponent,cell_data,trian,ds::DomainStyle)
  @notimplemented
end
function FESpaces.similar_fe_basis(
  f::MultiFieldFEBasisComponent,cell_data,trian,bs::BasisStyle,ds::DomainStyle)
  @notimplemented
end

for fun in (:gradient,:DIV,:∇∇)
  @eval begin
    function $fun(f::MultiFieldFEBasisComponent)
      g = $fun(f.single_field)
      MultiFieldFEBasisComponent(g,f.fieldid,f.nfields)
    end
  end
end

function CellData.change_domain(
  a::MultiFieldFEBasisComponent,
  tdomain::DomainStyle)
  sf = change_domain(a.single_field,tdomain)
  MultiFieldFEBasisComponent(sf,a.fieldid,a.nfields)
end

function CellData.change_domain(
  a::MultiFieldFEBasisComponent,
  ttrian::Triangulation,
  tdomain::DomainStyle)
  sf = change_domain(a.single_field,ttrian,tdomain)
  MultiFieldFEBasisComponent(sf,a.fieldid,a.nfields)
end

function MultiFieldFEBasisComponent(
  single_field::CellFieldAt{S,<:SingleFieldFEBasis},
  fieldid::Integer,
  nfields::Integer) where S

  sf = single_field.parent
  mf = MultiFieldFEBasisComponent(sf,fieldid,nfields)
  CellFieldAt{S}(mf)
end

function CellData.change_domain(
  a::CellFieldAt{S,<:MultiFieldFEBasisComponent},
  tdomain::DomainStyle) where S
  mf = a.parent
  sfin = CellFieldAt{S}(mf.single_field)
  sfout = change_domain(sfin,tdomain)
  MultiFieldFEBasisComponent(sfout,mf.fieldid,mf.nfields)
end

function CellData.change_domain(
  a::CellFieldAt{S,<:MultiFieldFEBasisComponent},
  ttrian::Triangulation,
  tdomain::DomainStyle) where S
  mf = a.parent
  sfin = CellFieldAt{S}(mf.single_field)
  sfout = change_domain(sfin,ttrian,tdomain)
  MultiFieldFEBasisComponent(sfout,mf.fieldid,mf.nfields)
end

function FESpaces.get_fe_basis(f::MultiFieldFESpace)
  nfields = length(f.spaces)
  all_febases = MultiFieldFEBasisComponent[]
  for field_i in 1:nfields
    dv_i = get_fe_basis(f.spaces[field_i])
    @assert BasisStyle(dv_i) == TestBasis()
    dv_i_b = MultiFieldFEBasisComponent(dv_i,field_i,nfields)
    push!(all_febases,dv_i_b)
  end
  MultiFieldCellField(all_febases)
end

function FESpaces.get_trial_fe_basis(f::MultiFieldFESpace)
  nfields = length(f.spaces)
  all_febases = MultiFieldFEBasisComponent[]
  for field_i in 1:nfields
    du_i = get_trial_fe_basis(f.spaces[field_i])
    @assert BasisStyle(du_i) == TrialBasis()
    du_i_b = MultiFieldFEBasisComponent(du_i,field_i,nfields)
    push!(all_febases,du_i_b)
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
    cell_values_field = lazy_map(a->a.array[i],cell_values)
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
  @check all(map(t->t===trian,trians)) msg
  get_cell_isconstrained(f,trian)
end

function FESpaces.get_cell_isconstrained(f::MultiFieldFESpace,trian::Triangulation)
  data = map(f.spaces) do space
    trian_i = get_triangulation(space)
    if is_change_possible(trian_i,trian)
      get_cell_isconstrained(space,trian)
    else
      Fill(false,num_cells(trian))
    end
  end
  lazy_map( (args...) -> +(args...)>0,  data...)
end

#function FESpaces.get_cell_isconstrained(f::MultiFieldFESpace,trian::SkeletonTriangulation)
#  plus = get_cell_isconstrained(f,trian.plus)
#  minus = get_cell_isconstrained(f,trian.minus)
#  lazy_map((l,r)-> l||r,plus,minus)
#end

function FESpaces.get_cell_is_dirichlet(f::MultiFieldFESpace)
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get
  the constrain flag for the corresponding cells.
  """
  trians = map(get_triangulation,f.spaces)
  trian = first(trians)
  @check all(map(t->t===trian,trians)) msg
  get_cell_is_dirichlet(f,trian)
end

function FESpaces.get_cell_is_dirichlet(f::MultiFieldFESpace,trian::Triangulation)
  data = map(f.spaces) do space
    trian_i = get_triangulation(space)
    if is_change_possible(trian_i,trian)
      get_cell_is_dirichlet(space,trian)
    else
      Fill(false,num_cells(trian))
    end
  end
  lazy_map( (args...) -> +(args...)>0,  data...)
end

#function FESpaces.get_cell_is_dirichlet(f::MultiFieldFESpace,trian::SkeletonTriangulation)
#  plus = get_cell_is_dirichlet(f,trian.plus)
#  minus = get_cell_is_dirichlet(f,trian.minus)
#  lazy_map((l,r)-> l||r,plus,minus)
#end

function FESpaces.get_cell_constraints(f::MultiFieldFESpace)
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get
  the constrains for the corresponding cells.
  """
  trians = map(get_triangulation,f.spaces)
  trian = first(trians)
  @check all(map(t->t===trian,trians)) msg
  get_cell_constraints(f,trian)
end

function FESpaces.get_cell_constraints(f::MultiFieldFESpace,trian::Triangulation)
  nfields = length(f.spaces)
  blockmask = [ is_change_possible(get_triangulation(Vi),trian) for Vi in f.spaces ]
  active_block_ids = findall(blockmask)
  active_block_data = Any[ get_cell_constraints(f.spaces[i],trian) for i in active_block_ids ]
  blockshape = (nfields,nfields)
  blockindices = [(i,i) for i in active_block_ids]
  lazy_map(BlockMap(blockshape,blockindices),active_block_data...)
end

#function FESpaces.get_cell_constraints(f::MultiFieldFESpace,trian::SkeletonTriangulation)
#  cell_values_plus = get_cell_constraints(f,trian.plus)
#  cell_values_minus = get_cell_constraints(f,trian.minus)
#  lazy_map(BlockMap((2,2),[(1,1),(2,2)]),cell_values_plus,cell_values_minus)
#end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace)
  msg = """\n
  This method does not make sense for multi-field
  since each field can be defined on a different triangulation.
  Pass a triangulation in the second argument to get the DOF ids
  on top of the corresponding cells.
  """
  trians = map(get_triangulation,f.spaces)
  trian = first(trians)
  @check all(map(t->t===trian,trians)) msg
  get_cell_dof_ids(f,trian)
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,trian::Triangulation)
  get_cell_dof_ids(f,trian,MultiFieldStyle(f))
end

#function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,trian::SkeletonTriangulation)
#  cell_values_plus = get_cell_dof_ids(f,trian.plus)
#  cell_values_minus = get_cell_dof_ids(f,trian.minus)
#  lazy_map(BlockMap(2,[1,2]),cell_values_plus,cell_values_minus)
#end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,::Triangulation,::MultiFieldStyle)
  @notimplemented
end

function FESpaces.get_cell_dof_ids(f::MultiFieldFESpace,trian::Triangulation,::ConsecutiveMultiFieldStyle)
  offsets = compute_field_offsets(f)
  nfields = length(f.spaces)
  blockmask = [ is_change_possible(get_triangulation(Vi),trian) for Vi in f.spaces ]
  active_block_ids = findall(blockmask)
  active_block_data = Any[]
  for i in active_block_ids
    cell_dofs_i = get_cell_dof_ids(f.spaces[i],trian)
    if i == 1
      push!(active_block_data,cell_dofs_i)
    else
      offset = Int32(offsets[i])
      o = Fill(offset,length(cell_dofs_i))
      cell_dofs_i_b = lazy_map(Broadcasting(_sum_if_first_positive),cell_dofs_i,o)
      push!(active_block_data,cell_dofs_i_b)
    end
  end
  lazy_map(BlockMap(nfields,active_block_ids),active_block_data...)
end

function _sum_if_first_positive(a,b)
  if a>0
    return a+b
  else
    return a
  end
end

function Arrays.return_value(k::Broadcasting{typeof(_sum_if_first_positive)},dofs::VectorBlock,o)
  evaluate(k,dofs,o)
end

function Arrays.return_cache(k::Broadcasting{typeof(_sum_if_first_positive)},dofs::VectorBlock,o)
  i = first(findall(dofs.touched))
  dofsi = dofs.array[i]
  ci = return_cache(k,dofsi,o)
  odofsi = evaluate!(ci,k,dofsi,o)
  array = Vector{typeof(odofsi)}(undef,length(dofs.array))
  c = Vector{typeof(ci)}(undef,length(dofs.array))
  for i in 1:length(dofs.array)
    if dofs.touched[i]
      c[i] = return_cache(k,dofs.array[i],o)
    end
  end
  ArrayBlock(array,dofs.touched), c
end

function Arrays.evaluate!(cache,k::Broadcasting{typeof(_sum_if_first_positive)},dofs::VectorBlock,o)
  r,c = cache
  for i in 1:length(dofs.array)
    if dofs.touched[i]
      r[i] = evaluate!(c[i],k,dofs.array[i],o)
    end
  end
  r
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

Base.length(m::MultiFieldFESpace) = length(m.spaces)

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
  interpolate!(objects,free_values,fe)
end

function FESpaces.interpolate!(objects,free_values::AbstractVector,fe::MultiFieldFESpace)
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
