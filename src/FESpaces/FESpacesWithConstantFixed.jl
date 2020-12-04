#TODO really needed ? It seems that we can return the given space if fix_constant == false
abstract type ConstantApproach end
struct FixConstant      <: ConstantApproach end
struct DoNotFixConstant <: ConstantApproach end

struct FESpaceWithConstantFixed{CA<:ConstantApproach,S<:SingleFieldFESpace} <: SingleFieldFESpace
  space::S
  dof_to_fix::Int
  @doc """
      FESpaceWithConstantFixed(space::SingleFieldFESpace, fix_constant::Bool,
      dof_to_fix::Int=num_free_dofs(space))
  """
  function FESpaceWithConstantFixed(
    space::SingleFieldFESpace,
    fix_constant::Bool,
    dof_to_fix::Int=num_free_dofs(space))

    S = typeof(space)
    if (fix_constant && num_dirichlet_dofs(space)==0)
      new{FixConstant,S}(space,dof_to_fix)
    else
      new{DoNotFixConstant,S}(space,dof_to_fix)
    end
  end
end

ConstantApproach(::FESpaceWithConstantFixed{CA}) where CA = CA()

# Genuine functions
function num_free_dofs(f::FESpaceWithConstantFixed{FixConstant})
  num_free_dofs(f.space)-1
end

function num_free_dofs(f::FESpaceWithConstantFixed{DoNotFixConstant})
  num_free_dofs(f.space)
end

function get_cell_dof_ids(f::FESpaceWithConstantFixed{FixConstant})
  cell_dofs = get_cell_dof_ids(f.space)
  CellDofIdsWithDofFixed(cell_dofs,f.dof_to_fix)
end

function get_cell_dof_ids(f::FESpaceWithConstantFixed{DoNotFixConstant})
  get_cell_dof_ids(f.space)
end

num_dirichlet_dofs(f::FESpaceWithConstantFixed{FixConstant}) = 1

num_dirichlet_dofs(f::FESpaceWithConstantFixed{DoNotFixConstant}) = 0

num_dirichlet_tags(f::FESpaceWithConstantFixed{FixConstant}) = 1

num_dirichlet_tags(f::FESpaceWithConstantFixed{DoNotFixConstant}) = 0

get_dirichlet_dof_tag(f::FESpaceWithConstantFixed{FixConstant}) = Int8[1,]

get_dirichlet_dof_tag(f::FESpaceWithConstantFixed{DoNotFixConstant}) = Int8[]

function scatter_free_and_dirichlet_values(f::FESpaceWithConstantFixed{FixConstant},fv,dv)
  @assert length(dv) == 1
  _dv = similar(dv,eltype(dv),0)
  _fv = VectorWithEntryInserted(fv,f.dof_to_fix,dv[1])
  scatter_free_and_dirichlet_values(f.space,_fv,_dv)
end

function scatter_free_and_dirichlet_values(f::FESpaceWithConstantFixed{DoNotFixConstant},fv,dv)
  @assert length(dv) == 0
 scatter_free_and_dirichlet_values(f.space,fv,dv)
end

function gather_free_and_dirichlet_values(f::FESpaceWithConstantFixed{FixConstant},cv)
  _fv, _dv = gather_free_and_dirichlet_values(f.space,cv)
  @assert length(_dv) == 0
  fv = VectorWithEntryRemoved(_fv,f.dof_to_fix)
  dv = _fv[f.dof_to_fix:f.dof_to_fix]
  (fv, dv)
end

function gather_free_and_dirichlet_values(f::FESpaceWithConstantFixed{DoNotFixConstant},cv)
  gather_free_and_dirichlet_values(f.space,cv)
end

function gather_free_and_dirichlet_values!(fv,dv,f::FESpaceWithConstantFixed{FixConstant},cv)
  _fv, _dv = gather_free_and_dirichlet_values(f.space,cv)
  @assert length(_dv) == 0
  fv    .= VectorWithEntryRemoved(_fv,f.dof_to_fix)
  dv[1]  = _fv[f.dof_to_fix]
  (fv, dv)
end

function gather_free_and_dirichlet_values!(fv,dv,f::FESpaceWithConstantFixed{DoNotFixConstant},cv)
  gather_free_and_dirichlet_values(f.space,cv)
end

function TrialFESpace(f::FESpaceWithConstantFixed{CA}) where CA
  U = TrialFESpace(f.space)
  FESpaceWithConstantFixed(U,CA==FixConstant,f.dof_to_fix)
end

# Delegated functions

function get_cell_shapefuns(f::FESpaceWithConstantFixed)
  get_cell_shapefuns(f.space)
end

function get_cell_shapefuns_trial(f::FESpaceWithConstantFixed)
  get_cell_shapefuns_trial(f.space)
end

function get_cell_dof_basis(f::FESpaceWithConstantFixed)
  get_cell_dof_basis(f.space)
end

get_vector_type(f::FESpaceWithConstantFixed) = get_vector_type(f.space)

get_dof_value_type(f::FESpaceWithConstantFixed) = get_dof_value_type(f.space)

CellField(t::FESpaceWithConstantFixed,cell_vals) = CellField(t.space,cell_vals)

ConstraintStyle(::Type{<:FESpaceWithConstantFixed{CA,S}}) where {CA,S} = ConstraintStyle(S)

get_triangulation(f::FESpaceWithConstantFixed) = get_triangulation(f.space)

# Helpers

struct CellDofIdsWithDofFixed{A<:AbstractArray} <: AbstractVector{Vector{Int32}}
  cell_dofs::A
  dof_to_fix::Int
end

Base.size(a::CellDofIdsWithDofFixed) = (length(a.cell_dofs),)

Base.IndexStyle(::Type{<:CellDofIdsWithDofFixed}) = IndexLinear()

function Base.getindex(a::CellDofIdsWithDofFixed,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function array_cache(a::CellDofIdsWithDofFixed)
  @assert eltype(a.cell_dofs) == Vector{Int32}
  b = testitem(a.cell_dofs)
  c = CachedArray(b)
  cache = array_cache(a.cell_dofs)
  (c, cache)
end

@inline function getindex!(d,a::CellDofIdsWithDofFixed,i::Integer)
  c, cache = d
  b = getindex!(cache,a.cell_dofs,i)
  setsize!(c,size(b))
  r = c.array
  for j in 1:length(b)
    bj = b[j]
    if bj == a.dof_to_fix
      rj = -1
    elseif bj < a.dof_to_fix
      rj = bj
    else
      rj = bj-1
    end
    r[j] = Int32(rj)
  end
  r
end
