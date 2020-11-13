
abstract type FEFunction <: CellField end

"""
"""
function get_free_values(f::FEFunction)
  @abstractmethod
end

"""
"""
function get_cell_dof_values(f::FEFunction)
  @abstractmethod
end

function get_cell_dof_values(f::FEFunction,cellids::AbstractArray)
  cell_values = get_cell_dof_values(f)
  lazy_map(Reindex(cell_values),cellids)
end

function get_cell_dof_values(f::FEFunction,cellids::SkeletonPair)
  cell_values_plus = get_cell_dof_values(f,cellids.plus)
  cell_values_minus = get_cell_dof_values(f,cellids.minus)
  cell_axes_plus = lazy_map(axes,cell_values_plus)
  cell_axes_minus = lazy_map(axes,cell_values_minus)
  cell_axes = lazy_map(cell_axes_plus,cell_axes_minus) do axp, axm
    (append_ranges([axp[1],axm[1]]),)
  end
  lazy_map(BlockArrayCooMap((2,),[(1,),(2,)]),cell_axes,cell_values_plus,cell_values_minus)
end

"""
"""
function get_fe_space(f::FEFunction)
  @abstractmethod
end

"""
"""
function test_fe_function(f::FEFunction)
  free_values = get_free_values(f)
  fe_space = get_fe_space(f)
  @test length(free_values) == num_free_dofs(fe_space)
  cell_values = get_cell_dof_values(f)
  trian = get_triangulation(f)
  @test length(cell_values) == num_cells(trian)
end

abstract type FESpace <: GridapType end

# Minimal FE interface (used by FEOperator)

"""
"""
function num_free_dofs(f::FESpace)
  @abstractmethod
end

"""
"""
function zero_free_values(f::FESpace)
  V = get_vector_type(f)
  allocate_vector(V,num_free_dofs(f))
end

function get_vector_type(fs::FESpace)
  @abstractmethod
end

"""
"""
function FEFunction(fe::FESpace, free_values)
  @abstractmethod
end

function CellField(fe::FESpace,cell_vals)
  @abstractmethod
end

function get_triangulation(fe::FESpace)
  @abstractmethod
end

# TODO this is quite hacky. Only needed by the zero mean space
function EvaluationFunction(fe::FESpace, free_values)
  FEFunction(fe,free_values)
end

"""
"""
function Base.zero(f::FESpace)
  free_values = zero_free_values(f)
  FEFunction(f,free_values)
end

function get_dof_value_type(f::FESpace)
  get_dof_value_type(get_cell_shapefuns(f),get_cell_dof_basis(f))
end

function get_dof_value_type(cell_shapefuns::CellField,cell_dof_basis::CellDof)
  cell_dof_values = cell_dof_basis(cell_shapefuns)
  eltype(eltype(cell_dof_values))
end

# Extended FEInterface used by Assemblers

"""
"""
function get_cell_dof_ids(f::FESpace)
  @abstractmethod
end

function get_cell_dof_ids(f::FESpace,cellids::AbstractArray)
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Reindex(cell_dof_ids),cellids)
end

function get_cell_dof_ids(f::FESpace,cellids::SkeletonPair)
  cell_ids_plus = get_cell_dof_ids(f,cellids.plus)
  cell_ids_minus = get_cell_dof_ids(f,cellids.minus)
  cell_axes_plus = lazy_map(axes,cell_ids_plus)
  cell_axes_minus = lazy_map(axes,cell_ids_minus)
  cell_axes = lazy_map(cell_axes_plus,cell_axes_minus) do axp, axm
    (append_ranges([axp[1],axm[1]]),)
  end
  lazy_map(BlockArrayCooMap((2,),[(1,),(2,)]),cell_axes,cell_ids_plus,cell_ids_minus)
end

"""
"""
function get_cell_shapefuns(f::FESpace)
  @abstractmethod
end

"""
"""
function get_cell_dof_basis(f::FESpace)
  @abstractmethod
end

function get_cell_shapefuns_trial(f::FESpace)
  v = get_cell_shapefuns(f)
  cell_v = get_cell_data(v)
  cell_u = lazy_map(transpose,cell_v)
  FEBasis(cell_u,get_triangulation(v),TrialBasis(),DomainStyle(v))
end

# Skeleton-related

abstract type BasisStyle end
struct TrialBasis <: BasisStyle end
struct TestBasis <: BasisStyle end

struct FEBasis{BS<:BasisStyle,DS<:DomainStyle} <: CellField
  cell_basis::AbstractArray{<:AbstractArray{<:Field}}
  trian::Triangulation
  basis_style::BS
  domain_style::DS
end

get_cell_data(f::FEBasis) = f.cell_basis
get_triangulation(f::FEBasis) = f.trian
BasisStyle(::Type{FEBasis{BS,DS}}) where {BS,DS} = BS()
BasisStyle(::T) where T <: FEBasis = BasisStyle(T)
DomainStyle(::Type{FEBasis{BS,DS}}) where {BS,DS} = DS()
function gradient(a::FEBasis)
  f = GenericCellField(a.cell_basis,a.trian,a.domain_style)
  FEBasis(get_cell_data(gradient(f)),a.trian,a.basis_style,a.domain_style)
end
function ∇∇(a::FEBasis)
  f = GenericCellField(a.cell_basis,a.trian,a.domain_style)
  FEBasis(get_cell_data(∇∇(f)),a.trian,a.basis_style,a.domain_style)
end

function change_domain_skeleton(a::FEBasis,trian::SkeletonTriangulation,target_domain::DomainStyle)
  a_on_plus_trian = change_domain(a,trian.plus,target_domain)
  a_on_minus_trian = change_domain(a,trian.minus,target_domain)
  pair_in = SkeletonPair(get_cell_data(a_on_plus_trian),get_cell_data(a_on_minus_trian))
  pair_out = _fix_cell_basis_dofs_at_skeleton(pair_in,BasisStyle(a))
  plus = GenericCellField(pair_out.plus,trian,target_domain)
  minus = GenericCellField(pair_out.minus,trian,target_domain)
  plus, minus
end

function _fix_cell_basis_dofs_at_skeleton(pair,::TestBasis)
  cell_axes_plus = lazy_map(axes,pair.plus)
  cell_axes_minus = lazy_map(axes,pair.minus)
  function skeleton_axes_test(axp, axm)
    (append_ranges([axp[1],axm[1]]),)
  end
  cell_axes = lazy_map(skeleton_axes_test,cell_axes_plus,cell_axes_minus)
  plus = lazy_map(BlockFieldArrayCooMap((2,),[(1,)]),cell_axes,pair.plus)
  minus = lazy_map(BlockFieldArrayCooMap((2,),[(2,)]),cell_axes,pair.minus)
  SkeletonPair(plus,minus)
end

function _fix_cell_basis_dofs_at_skeleton(pair,::TrialBasis)
  cell_axes_plus = lazy_map(axes,pair.plus)
  cell_axes_minus = lazy_map(axes,pair.minus)
  function skeleton_axes_trial(axp, axm)
    r2 = append_ranges([axp[2],axm[2]])
    r1 = similar_range(r2,1)
    (r1,r2)
  end
  cell_axes = lazy_map(skeleton_axes_trial,cell_axes_plus,cell_axes_minus)
  plus = lazy_map(BlockFieldArrayCooMap((1,2),[(1,1)]),cell_axes,pair.plus)
  minus = lazy_map(BlockFieldArrayCooMap((1,2),[(1,2)]),cell_axes,pair.minus)
  SkeletonPair(plus,minus)
end

# Constraint-related

abstract type ConstraintStyle end
struct Constrained <: ConstraintStyle end
struct UnConstrained <: ConstraintStyle end

ConstraintStyle(::Type{<:FESpace}) = @abstractmethod
ConstraintStyle(::T) where T<:FESpace = ConstraintStyle(T)

has_constraints(::Type{T}) where T <:FESpace = ConstraintStyle(T) == Constrained()
has_constraints(::T) where T <: FESpace = has_constraints(T)

function get_cell_constraints(f::FESpace)
  get_cell_constraints(f,ConstraintStyle(f))
end

function get_cell_constraints(f,::UnConstrained)
  cell_axes = lazy_map(axes,get_cell_data(get_cell_shapefuns(f)))
  identity_constraints(cell_axes)
end

function get_cell_constraints(f,::Constrained)
  @abstractmethod
end

function get_cell_constraints(f::FESpace,cellids::AbstractArray)
  lazy_map(Reindex(get_cell_constraints(f)),cellids)
end

function get_cell_constraints(f::FESpace,cellids::SkeletonPair)
  cell_constraints_plus = get_cell_constraints(f,cellids.plus)
  cell_constraints_minus = get_cell_constraints(f,cellids.minus)
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

function get_cell_isconstrained(f::FESpace)
  get_cell_isconstrained(f,ConstraintStyle(f))
end

function get_cell_isconstrained(f,::UnConstrained)
  Fill(false,length(get_cell_dof_ids(f)))
end

function get_cell_isconstrained(f,::Constrained)
  @abstractmethod
end

function get_cell_isconstrained(f::FESpace,cellids::AbstractArray)
  lazy_map(Reindex(get_cell_isconstrained(f)),cellids)
end

function get_cell_isconstrained(f::FESpace,cellids::SkeletonPair)
  plus = get_cell_isconstrained(f,cellids.plus)
  minus = get_cell_isconstrained(f,cellids.minus)
  lazy_map((l,r)-> l||r,plus,minus)
end

function attach_constraints_rows(f::FESpace,cellarr,cellids)
  attach_constraints_rows(f,cellarr,cellids,ConstraintStyle(f))
end

function attach_constraints_rows(f::FESpace,cellarr,cellids,::UnConstrained)
  cellarr
end

function attach_constraints_rows(f::FESpace,cellarr,cellids,::Constrained)
  cellconstr = get_cell_constraints(f,cellids)
  cellmask = get_cell_isconstrained(f,cellids)
  attach_constraints_rows(cellarr,cellconstr,cellmask)
end

function attach_constraints_cols(f::FESpace,cellarr,cellids)
  attach_constraints_cols(f,cellarr,cellids,ConstraintStyle(f))
end

function attach_constraints_cols(f::FESpace,cellarr,cellids,::UnConstrained)
  cellarr
end

function attach_constraints_cols(f::FESpace,cellarr,cellids,::Constrained)
  cellconstr = get_cell_constraints(f,cellids)
  cellmask = get_cell_isconstrained(f,cellids)
  attach_constraints_cols(cellarr,cellconstr,cellmask)
end

"""
"""
function test_fe_space(f::FESpace)
  trian = get_triangulation(f)
  @test isa(trian,Triangulation)
  free_values = zero_free_values(f)
  @test length(free_values) == num_free_dofs(f)
  V = get_vector_type(f)
  @test typeof(free_values) == V
  fe_function = FEFunction(f,free_values)
  test_fe_function(fe_function)
  fe_basis = get_cell_shapefuns(f)
  @test isa(has_constraints(f),Bool)
  @test isa(has_constraints(typeof(f)),Bool)
  @test length(get_cell_dof_ids(f)) == num_cells(fe_basis)
  @test length(get_cell_constraints(f)) == num_cells(fe_basis)
  @test length(get_cell_isconstrained(f)) == num_cells(fe_basis)
  @test CellField(f,get_cell_dof_ids(f)) != nothing
end

function test_fe_space(f::FESpace,matvecdata,matdata,vecdata)
  test_fe_space(f)

  cellmat, cellidsrows, cellidscols = matdata
  cm = attach_constraints_cols(f,cellmat,cellidscols)
  if ! has_constraints(f)
    @test cm === cellmat
  end
  cm = attach_constraints_rows(f,cellmat,cellidsrows)
  if ! has_constraints(f)
    @test cm === cellmat
  end

  cellvec, cellidsrows = vecdata
  cv = attach_constraints_rows(f,cellvec,cellidsrows)
  if ! has_constraints(f)
    @test cv === cellvec
  end

  cellmatvec, cellidsrows, cellidscols = matvecdata
  cmv = attach_constraints_cols(f,cellmatvec,cellidscols)
  if ! has_constraints(f)
    @test cmv === cellmatvec
  end
  cmv = attach_constraints_rows(f,cellmatvec,cellidsrows)
  if ! has_constraints(f)
    @test cmv === cellmatvec
  end

end

