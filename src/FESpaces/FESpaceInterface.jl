
abstract type FEFunction <: CellField end

"""
"""
function get_free_dof_values(f::FEFunction)
  @abstractmethod
end

"""
"""
function get_cell_dof_values(f::FEFunction)
  @abstractmethod
end

function get_cell_dof_values(f::FEFunction,ttrian::Triangulation)
  get_cell_fe_data(get_cell_dof_values,f,ttrian)
end

function get_cell_fe_data(fun,f,ttrian)
  sface_to_data = fun(f)
  strian = get_triangulation(f)
  if strian === ttrian
    return sface_to_data
  end
  @check is_change_possible(strian,ttrian)
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  get_cell_fe_data(fun,sface_to_data,sglue,tglue)
end

function get_cell_fe_data(fun,sface_to_data,sglue::FaceToFaceGlue,tglue::FaceToFaceGlue)
  mface_to_sface = sglue.mface_to_tface
  tface_to_mface = tglue.tface_to_mface
  mface_to_data = extend(sface_to_data,mface_to_sface)
  tface_to_data = lazy_map(Reindex(mface_to_data),tface_to_mface)
  tface_to_data
end

function get_cell_fe_data(fun,sface_to_data,sglue::FaceToFaceGlue,tglue::SkeletonPair)
  plus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.plus)
  minus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.minus)
  lazy_map(BlockMap(2,[1,2]),plus,minus)
end

"""
"""
function get_fe_space(f::FEFunction)
  @abstractmethod
end

function get_cell_is_dirichlet(f::FEFunction)
  get_cell_is_dirichlet(get_fe_space(f))
end

function get_cell_is_dirichlet(f::FEFunction,ttrian::Triangulation)
  get_cell_is_dirichlet(get_fe_space(f),ttrian)
end

"""
"""
function test_fe_function(f::FEFunction)
  trian = get_triangulation(f)
  free_values = get_free_dof_values(f)
  fe_space = get_fe_space(f)
  @test length(free_values) == num_free_dofs(fe_space)
  cell_values = get_cell_dof_values(f,trian)
  @test length(cell_values) == num_cells(trian)
end

"""
    abstract type FESpace <: GridapType end

Abstract finite element space.
"""
abstract type FESpace <: GridapType end

# Minimal FE interface (used by FEOperator)

"""
    get_free_dof_ids(f::FESpace)
"""
function get_free_dof_ids(f::FESpace)
  @abstractmethod
end

"""
    num_free_dofs(f::FESpace) = length(get_free_dof_ids(f))
"""
num_free_dofs(f::FESpace) = length(get_free_dof_ids(f))

"""
    zero_free_values(f::FESpace)
"""
function zero_free_values(f::FESpace)
  V = get_vector_type(f)
  free_values = allocate_vector(V,get_free_dof_ids(f))
  fill!(free_values,zero(eltype(V)))
  return free_values
end

"""
    get_vector_type(fs::FESpace)
"""
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
"""
    EvaluationFunction(fe::FESpace, free_values) = FEFunction(fe,free_values)
"""
function EvaluationFunction(fe::FESpace, free_values)
  FEFunction(fe,free_values)
end

"""
"""
function Base.zero(f::FESpace)
  free_values = zero_free_values(f)
  EvaluationFunction(f,free_values)
end

function get_dof_value_type(f::FESpace)
  get_dof_value_type(get_fe_basis(f),get_fe_dof_basis(f))
end

"""
"""
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

function get_cell_dof_ids(f::FESpace,ttrian::Triangulation)
  get_cell_fe_data(get_cell_dof_ids,f,ttrian)
end

"""
"""
function get_fe_basis(f::FESpace)
  @abstractmethod
end

function get_cell_shapefuns(f::FESpace)
  msg = "get_cell_shapefuns has been removed, use get_fe_basis instead"
  error(msg)
end

"""
"""
function get_fe_dof_basis(f::FESpace)
  @abstractmethod
end

function get_cell_dof_basis(f::FESpace)
  msg = "get_cell_dof_basis has been removed, use get_fe_dof_basis instead"
  error(msg)
end

"""
"""
function get_trial_fe_basis(f::FESpace)
  v = get_fe_basis(f)
  cell_v = get_data(v)
  cell_u = lazy_map(transpose,cell_v)
  similar_fe_basis(v,cell_u,get_triangulation(v),TrialBasis(),DomainStyle(v))
end

function get_cell_shapefuns_trial(f::FESpace)
  msg = "get_cell_shapefuns_trial has been removed, use get_trial_fe_basis instead"
  error(msg)
end

# Basis related

"""
    abstract type BasisStyle

Trait for trial or test [`FEBasis`](@ref), the subtypes are the structs `TrialBasis` and `TestBasis`.
"""
abstract type BasisStyle end
struct TrialBasis <: BasisStyle end
struct TestBasis <: BasisStyle end

"""
    abstract type FEBasis <: CellField

Has traits [BasisStyle](@ref) and [DomainStyle](@ref).
"""
abstract type FEBasis <: CellField end
BasisStyle(::T) where T <: FEBasis = BasisStyle(T)
DomainStyle(::T) where T <: FEBasis = DomainStyle(T)

struct SingleFieldFEBasis{BS<:BasisStyle,DS<:DomainStyle} <: FEBasis
  cell_basis::AbstractArray
  trian::Triangulation
  basis_style::BS
  domain_style::DS
  function SingleFieldFEBasis(
    cell_basis::AbstractArray,
    trian::Triangulation,
    basis_style::BasisStyle,
    domain_style::DomainStyle)

    BS = typeof(basis_style)
    DS = typeof(domain_style)
    new{BS,DS}(Fields.MemoArray(cell_basis),trian,basis_style,domain_style)
  end
end

get_data(f::SingleFieldFEBasis) = f.cell_basis
get_triangulation(f::SingleFieldFEBasis) = f.trian
BasisStyle(::Type{SingleFieldFEBasis{BS,DS}}) where {BS,DS} = BS()
DomainStyle(::Type{SingleFieldFEBasis{BS,DS}}) where {BS,DS} = DS()
function CellData.similar_cell_field(f::SingleFieldFEBasis,cell_data,trian,ds::DomainStyle)
  SingleFieldFEBasis(cell_data,trian,BasisStyle(f),ds)
end
function similar_fe_basis(f::SingleFieldFEBasis,cell_data,trian,bs::BasisStyle,ds::DomainStyle)
  SingleFieldFEBasis(cell_data,trian,bs,ds)
end

for fun in (:change_domain_ref_ref,:change_domain_phys_phys)
  @eval begin

    function $fun(
      a::CellData.CellFieldAt{S,<:FEBasis} where S,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::SkeletonPair)
      a_on_plus_trian = $fun(a.parent,ttrian,sglue,tglue.plus)
      a_on_minus_trian = $fun(a.parent,ttrian,sglue,tglue.minus)
      pair_in = SkeletonPair(get_data(a_on_plus_trian),get_data(a_on_minus_trian))
      pair_out = _fix_cell_basis_dofs_at_skeleton(pair_in,BasisStyle(a.parent))
      if isa(a,CellData.CellFieldAt{:plus})
        return CellData.similar_cell_field(a.parent,pair_out.plus,ttrian,DomainStyle(a.parent))
      elseif isa(a,CellData.CellFieldAt{:minus})
        return CellData.similar_cell_field(a.parent,pair_out.minus,ttrian,DomainStyle(a.parent))
      else
        @unreachable
      end
    end

  end
end

function _fix_cell_basis_dofs_at_skeleton(pair,::TestBasis)
  plus = lazy_map(BlockMap(2,1),pair.plus)
  minus = lazy_map(BlockMap(2,2),pair.minus)
  SkeletonPair(plus,minus)
end

function _fix_cell_basis_dofs_at_skeleton(pair,::TrialBasis)
  plus = lazy_map(BlockMap((1,2),1),pair.plus)
  minus = lazy_map(BlockMap((1,2),2),pair.minus)
  SkeletonPair(plus,minus)
end

# Constraint-related

"""
    abstract type ConstraintStyle

Trait for (un)constrained [`FESpace`](@ref)s, the subtypes are the structs
[`Constrained`](@ref) and [`UnConstrained`](@ref).
"""
abstract type ConstraintStyle end
"""
    struct Constrained <: ConstraintStyle
"""
struct Constrained <: ConstraintStyle end
"""
    struct UnConstrained <: ConstraintStyle
"""
struct UnConstrained <: ConstraintStyle end

ConstraintStyle(::Type{<:FESpace}) = @abstractmethod
ConstraintStyle(::T) where T<:FESpace = ConstraintStyle(T)

"""
    has_constraints(::Type{<FESpace})
    has_constraints(::FESpace)

Return true if the `FESpace` (type) is [`Constrained`](@ref).
"""
has_constraints(::Type{T}) where T <:FESpace = ConstraintStyle(T) == Constrained()
has_constraints(::T) where T <: FESpace = has_constraints(T)

"""
"""
function get_cell_constraints(f::FESpace)
  get_cell_constraints(f,ConstraintStyle(f))
end

function get_cell_constraints(f::FESpace,::UnConstrained)
  cell_ids = get_cell_dof_ids(f)
  cell_axes = lazy_map(axes,cell_ids)
  identity_constraints(cell_axes)
end

function get_cell_constraints(f::FESpace,::Constrained)
  @abstractmethod
end

function get_cell_constraints(f::FESpace,ttrian::Triangulation)
  get_cell_fe_data(get_cell_constraints,f,ttrian)
end

function get_cell_fe_data(fun::typeof(get_cell_constraints),sface_to_data,sglue::FaceToFaceGlue,tglue::SkeletonPair)
  plus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.plus)
  minus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.minus)
  lazy_map(BlockMap((2,2),[(1,1),(2,2)]),plus,minus)
end

"""
"""
function get_cell_isconstrained(f::FESpace)
  get_cell_isconstrained(f,ConstraintStyle(f))
end

function get_cell_isconstrained(f::FESpace,::UnConstrained)
  trian = get_triangulation(f)
  Fill(false,num_cells(trian))
end

function get_cell_isconstrained(f::FESpace,::Constrained)
  @abstractmethod
end

function get_cell_isconstrained(f::FESpace,ttrian::Triangulation)
  get_cell_fe_data(get_cell_isconstrained,f,ttrian)
end

function get_cell_fe_data(
  fun::typeof(get_cell_isconstrained),sface_to_data,sglue::FaceToFaceGlue,tglue::SkeletonPair)
  plus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.plus)
  minus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.minus)
  lazy_map((l,r)-> l||r,plus,minus)
end

function attach_constraints_rows(f::FESpace,cellarr,ttrian::Triangulation)
  attach_constraints_rows(f,cellarr,ttrian,ConstraintStyle(f))
end

function attach_constraints_rows(f::FESpace,cellarr,ttrian,::UnConstrained)
  cellarr
end

function attach_constraints_rows(f::FESpace,cellarr,ttrian,::Constrained)
  cellconstr = get_cell_constraints(f,ttrian)
  cellmask = get_cell_isconstrained(f,ttrian)
  attach_constraints_rows(cellarr,cellconstr,cellmask)
end

function attach_constraints_cols(f::FESpace,cellarr,ttrian::Triangulation)
  attach_constraints_cols(f,cellarr,ttrian,ConstraintStyle(f))
end

function attach_constraints_cols(f::FESpace,cellarr,ttrian,::UnConstrained)
  cellarr
end

function attach_constraints_cols(f::FESpace,cellarr,ttrian,::Constrained)
  cellconstr = get_cell_constraints(f,ttrian)
  cellmask = get_cell_isconstrained(f,ttrian)
  attach_constraints_cols(cellarr,cellconstr,cellmask)
end

"""
"""
function get_cell_is_dirichlet(f::FESpace)
  trian = get_triangulation(f)
  Fill(true,num_cells(trian))
end

function get_cell_is_dirichlet(f::FESpace,ttrian::Triangulation)
  get_cell_fe_data(get_cell_is_dirichlet,f,ttrian)
end

function get_cell_fe_data(
  fun::typeof(get_cell_is_dirichlet),sface_to_data,sglue::FaceToFaceGlue,tglue::SkeletonPair)
  plus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.plus)
  minus = get_cell_fe_data(fun,sface_to_data,sglue,tglue.minus)
  lazy_map((l,r)-> l||r,plus,minus)
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
  fe_basis = get_fe_basis(f)
  @test isa(has_constraints(f),Bool)
  @test isa(has_constraints(typeof(f)),Bool)
  @test length(get_cell_dof_ids(f,trian)) == num_cells(fe_basis)
  @test length(get_cell_constraints(f,trian)) == num_cells(fe_basis)
  @test length(get_cell_isconstrained(f,trian)) == num_cells(fe_basis)
  @test CellField(f,get_cell_dof_ids(f,trian)) != nothing
end

function test_fe_space(f::FESpace,cell_matvec,cell_mat,cell_vec,trian)
  test_fe_space(f)

  cm = attach_constraints_cols(f,cell_mat,trian)
  if ! has_constraints(f)
    @test cm === cell_mat
  end
  cm = attach_constraints_rows(f,cell_mat,trian)
  if ! has_constraints(f)
    @test cm === cell_mat
  end

  cv = attach_constraints_rows(f,cell_vec,trian)
  if ! has_constraints(f)
    @test cv === cell_vec
  end

  cmv = attach_constraints_cols(f,cell_matvec,trian)
  if ! has_constraints(f)
    @test cmv === cell_matvec
  end
  cmv = attach_constraints_rows(f,cell_matvec,trian)
  if ! has_constraints(f)
    @test cmv === cell_matvec
  end

end
