
struct ExtendedFEFunctionValue{T,A,B} <: AbstractVector{Vector{T}}
  cell_to_val::A
  oldcell_to_cell::Vector{Int}
  oldcell_to_points::B
end

Base.size(a::ExtendedFEFunctionValue) = (length(a.oldcell_to_cell),)

Base.IndexStyle(::Type{<:ExtendedFEFunctionValue}) = IndexLinear()

function Base.getindex(a::ExtendedFEFunctionValue,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function array_cache(a::ExtendedFEFunctionValue)
  c1 = array_cache(a.oldcell_to_cell)
  c2 = array_cache(a.oldcell_to_points)
  val = testitem(a.cell_to_val)
  v = copy(val)
  ca = CachedArray(v)
  (ca,c1,c2)
end

@inline function getindex!(cache,a::ExtendedFEFunctionValue,oldcell)
  ca, c1, c2 = cache
  cell = a.oldcell_to_cell[oldcell]
  if cell == UNSET
    x = getindex!(c2,a.oldcell_to_points,oldcell)
    Q = length(x)
    setsize!(ca,(Q,))
    v = ca.array
    fill!(v,zero(eltype(v)))
    return v
  else
    return getindex!(c1,a.cell_to_val,cell)
  end
end

struct ExtendedFEFunctionField <: AbstractVector{UnimplementedField}
  cell_to_field
  cell_to_oldcell::Vector{Int}
  oldcell_to_cell::Vector{Int}
end

Base.size(a::ExtendedFEFunctionField) = (length(a.oldcell_to_cell),)

Base.IndexStyle(::Type{<:ExtendedFEFunctionField}) = IndexLinear()

function evaluate_field_array(a::ExtendedFEFunctionField,x::AbstractArray)
  _x = reindex(x,a.cell_to_oldcell)
  cell_to_val = evaluate_field_array(a.cell_to_field,_x)
  ExtendedFEFunctionValue(cell_to_val,oldcell_to_cell,x)
end

# TODO This is temporary
using Gridap.Arrays: Reindexed

"""
"""
struct ExtendedFESpace <: SingleFieldFESpace
  space::SingleFieldFESpace
  trian::RestrictedTriangulation
end

constraint_style(::Type{ExtendedFESpace}) = Val{false}()

function get_cell_dofs(f::ExtendedFESpace)
  Reindexed(get_cell_dofs(f.space),f.trian.oldcell_to_cell)
end

function get_cell_basis(f::ExtendedFESpace)
  cell_basis = get_cell_basis(f.space)
  a = get_array(cell_basis)
  array = Reindexed(a,f.trian.oldcell_to_cell)
  b = get_cell_map(cell_basis)
  cm = get_cell_map(f.trian.oldtrian)
  trial_style = TrialStyle(cell_basis)
  GenericCellBasis(trial_style,array,cm)
end

function get_cell_dof_basis(f::ExtendedFESpace)
  Reindexed(get_cell_dof_basis(f.space),f.trian.oldcell_to_cell)
end

function FEFunction(
  f::ExtendedFESpace, free_values::AbstractVector, dirichlet_values::AbstractVector)

  cell_uh = FEFunction(f.space,free_values,dirichlet_values)
  cell_field = ExtendedFEFunctionField(
    get_array(cell_uh),f.trian.cell_to_oldcell,f.trian.oldcell_to_cell)
  cell_vals = scatter_free_and_dirichlet_values(f,free_values,dirichlet_values)
  SingleFieldFEFunction(cell_field,cell_vals,free_values,dirichlet_values,f)
end

function scatter_free_and_dirichlet_values(f::ExtendedFESpace,fv,dv)
  a = scatter_free_and_dirichlet_values(f.space,fv,dv)
  Reindexed(a,f.trian.oldcell_to_cell)
end

function gather_free_and_dirichlet_values(f::ExtendedFESpace,cv)
  _cv = reindex(cv,f.trian.cell_to_oldcell)
  gather_free_and_dirichlet_values(f.space,_cv)
end





function num_free_dofs(f::ExtendedFESpace)
  num_free_dofs(f.space)
end

function zero_free_values(::Type{T},f::ExtendedFESpace) where T
  zeros(T,num_free_dofs(f))
end

function num_dirichlet_dofs(f::ExtendedFESpace)
  num_dirichlet_dofs(f.space)
end

function zero_dirichlet_values(f::ExtendedFESpace)
  zero_dirichlet_values(f.space)
end

function num_dirichlet_tags(f::ExtendedFESpace)
  num_dirichlet_tags(f.space)
end

function get_dirichlet_dof_tag(f::ExtendedFESpace)
  get_dirichlet_dof_tag(f.space)
end

function TrialFESpace(f::ExtendedFESpace)
  U = TrialFESpace(f.space)
  ExtendedFESpace(U,f.trian)
end

