

struct ExtendedVector{T,A,B} <: AbstractVector{T}
  void_to_val::A
  cell_to_val::B
  oldcell_to_cell_or_void::Vector{Int}
  void_to_oldcell::Vector{Int}
  cell_to_oldcell::Vector{Int}
  function ExtendedVector(
    void_to_val::AbstractArray,
    cell_to_val::AbstractArray,
    oldcell_to_cell_or_void::Vector{Int},
    void_to_oldcell::Vector{Int},
    cell_to_oldcell::Vector{Int})

    A = typeof(void_to_val)
    B = typeof(cell_to_val)
    ai = testitem(void_to_val)
    bi = testitem(cell_to_val)
    T = eltype([ai,bi])

    new{T,A,B}(
      void_to_val,
      cell_to_val,
      oldcell_to_cell_or_void,
      void_to_oldcell,
      cell_to_oldcell)

  end
end

Base.size(a::ExtendedVector) = (length(a.oldcell_to_cell_or_void),)

Base.IndexStyle(::Type{<:ExtendedVector}) = IndexLinear()

function Base.getindex(a::ExtendedVector,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function array_cache(a::ExtendedVector)
  cv = array_cache(a.void_to_val)
  cc = array_cache(a.cell_to_val)
  (cv,cc)
end

@inline function getindex!(cache,a::ExtendedVector,oldcell)
  cv, cc = cache
  i = a.oldcell_to_cell_or_void[oldcell]
  if i > 0
    cell = i
    return getindex!(cc,a.cell_to_val,cell)
  else
    void = -i
    return getindex!(cv,a.void_to_val,void)
  end
end

function evaluate_field_array(a::ExtendedVector,x::AbstractArray)
  x_void = reindex(x,a.void_to_oldcell)
  x_cell = reindex(x,a.cell_to_oldcell)
  r_void = evaluate_field_array(a.void_to_val,x_void)
  r_cell = evaluate_field_array(a.cell_to_val,x_cell)
  ExtendedVector(
      r_void,
      r_cell,
      oldcell_to_cell_or_void,
      void_to_oldcell,
      cell_to_oldcell)
end

function reindex(a::ExtendedVector,ptrs)
  if a.cell_to_oldcell === ptrs || a.cell_to_oldcell == ptrs
    return a.cell_to_val
  elseif a.void_to_oldcell === ptrs || a.void_to_oldcell == ptrs
    return a.void_to_val
  else
    return Reindexed(a,ptrs)
  end
end

struct VoidBasis{T} <: Field end

function field_cache(f::VoidBasis{T},x) where T
  Q = length(x)
  v = zeros(T,(Q,0))
  CachedArray(v)
end

@inline function evaluate_field!(cache,f::VoidBasis{T},x) where T
  Q = length(x)
  setsize!(cache,(Q,0))
  cache.array
end

function field_gradient(f::VoidBasis{T}) where T
  VoidBasis{eltype(T)}()
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

