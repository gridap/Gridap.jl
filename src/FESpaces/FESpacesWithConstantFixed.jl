abstract type ConstantApproach end;
struct FixConstant      <: ConstantApproach end;
struct DoNotFixConstant <: ConstantApproach end;

#"""
#    struct FESpaceWithConstantFixed{CS,<:ConstantApproach} <: SingleFieldFESpace
#      space::SingleFieldFESpace
#      constraint_style::Val{CS}
#      dof_to_fix::Int
#    end
#"""
struct FESpaceWithConstantFixed{CS,CA<:ConstantApproach} <: SingleFieldFESpace
  space::SingleFieldFESpace
  constraint_style::Val{CS}
  dof_to_fix::Int
  @doc """
      FESpaceWithConstantFixed(space::SingleFieldFESpace, fix_constant::Bool,
      dof_to_fix::Int=num_free_dofs(space))
  """
  function FESpaceWithConstantFixed(space::SingleFieldFESpace,
     fix_constant::Bool,
     dof_to_fix::Int=num_free_dofs(space))
    cs = constraint_style(space)
    CS = get_val_parameter(cs)
    if (fix_constant && num_dirichlet_dofs(space)==0)
      new{CS,FixConstant}(space,cs,dof_to_fix)
    else
      new{CS,DoNotFixConstant}(space,cs,dof_to_fix)
    end
  end
end

const FESpaceWithLastDofRemoved{CS} = FESpaceWithConstantFixed{CS,FixConstant}
@deprecate FESpaceWithLastDofRemoved(space::SingleFieldFESpace) = FESpaceWithConstantFixed(space,true)

# Genuine functions
function num_free_dofs(f::FESpaceWithConstantFixed{CS,FixConstant}) where {CS}
  num_free_dofs(f.space)-1
end

function num_free_dofs(f::FESpaceWithConstantFixed{CS,DoNotFixConstant}) where {CS}
  num_free_dofs(f.space)
end

function zero_free_values(f::FESpaceWithConstantFixed)
  zeros(num_free_dofs(f))
end

function get_cell_dofs(f::FESpaceWithConstantFixed{CS,FixConstant}) where {CS}
  cell_dofs = get_cell_dofs(f.space)
  CellDofsWithDofFixed(cell_dofs,f.dof_to_fix)
end

function get_cell_dofs(f::FESpaceWithConstantFixed{CS,DoNotFixConstant}) where {CS}
  get_cell_dofs(f.space)
end


num_dirichlet_dofs(f::FESpaceWithConstantFixed{CS,FixConstant}) where {CS} = 1

num_dirichlet_dofs(f::FESpaceWithConstantFixed{CS,DoNotFixConstant}) where {CS} = 0

function zero_dirichlet_values(f::FESpaceWithConstantFixed)
  T = Float64 # TODO
  zeros(T,num_dirichlet_dofs(f))
end


num_dirichlet_tags(f::FESpaceWithConstantFixed{CS,FixConstant}) where {CS} = 1

num_dirichlet_tags(f::FESpaceWithConstantFixed{CS,DoNotFixConstant}) where {CS} = 0

get_dirichlet_dof_tag(f::FESpaceWithConstantFixed{CS,FixConstant}) where {CS} = Int8[1,]

get_dirichlet_dof_tag(f::FESpaceWithConstantFixed{CS,DoNotFixConstant}) where {CS} = Int8[]

function scatter_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,FixConstant},fv,dv) where {CS}
  @assert length(dv) == 1
  _dv = similar(dv,eltype(dv),0)
  _fv = VectorWithEntryInserted(fv,f.dof_to_fix,dv[1])
  scatter_free_and_dirichlet_values(f.space,_fv,_dv)
end

function scatter_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,DoNotFixConstant},fv,dv) where {CS}
  @assert length(dv) == 0
 scatter_free_and_dirichlet_values(f.space,fv,dv)
end

function gather_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,FixConstant},cv) where {CS}
  _fv, _dv = gather_free_and_dirichlet_values(f.space,cv)
  @assert length(_dv) == 0
  fv = VectorWithEntryRemoved(_fv,f.dof_to_fix)
  dv = _fv[f.dof_to_fix:f.dof_to_fix]
  (fv, dv)
end

function gather_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,DoNotFixConstant},cv) where {CS}
  gather_free_and_dirichlet_values(f.space,cv)
end

function gather_free_and_dirichlet_values!(
  fv,dv,f::FESpaceWithConstantFixed{CS,FixConstant},cv) where {CS}
  _fv, _dv = gather_free_and_dirichlet_values(f.space,cv)
  @assert length(_dv) == 0
  fv    .= VectorWithEntryRemoved(_fv,f.dof_to_fix)
  dv[1]  = _fv[f.dof_to_fix]
  (fv, dv)
end

function gather_free_and_dirichlet_values!(
  fv,dv,f::FESpaceWithConstantFixed{CS,DoNotFixConstant},cv) where {CS}
  gather_free_and_dirichlet_values(f.space,cv)
end

function TrialFESpace(f::FESpaceWithConstantFixed{CS,CA}) where {CS,CA}
  U = TrialFESpace(f.space)
  FESpaceWithConstantFixed(U,CA==FixConstant,f.dof_to_fix)
end

# Delegated functions

function get_cell_basis(f::FESpaceWithConstantFixed)
  get_cell_basis(f.space)
end

function get_cell_dof_basis(f::FESpaceWithConstantFixed)
  get_cell_dof_basis(f.space)
end

get_cell_axes(t::FESpaceWithConstantFixed)= get_cell_axes(t.space)

get_cell_axes_with_constraints(t::FESpaceWithConstantFixed)= get_cell_axes_with_constraints(t.space)

CellData.CellField(t::FESpaceWithConstantFixed,cell_vals) = CellField(t.space,cell_vals)

constraint_style(::Type{<:FESpaceWithConstantFixed{CS}}) where CS = Val{CS}()

# Helpers

struct CellDofsWithDofFixed{A<:AbstractArray} <: AbstractVector{Vector{Int}}
  cell_dofs::A
  dof_to_fix::Int
end

Base.size(a::CellDofsWithDofFixed) = (length(a.cell_dofs),)

Base.IndexStyle(::Type{<:CellDofsWithDofFixed}) = IndexLinear()

function Base.getindex(a::CellDofsWithDofFixed,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function array_cache(a::CellDofsWithDofFixed)
  @assert eltype(a.cell_dofs) == Vector{Int}
  b = testitem(a.cell_dofs)
  c = CachedArray(b)
  cache = array_cache(a.cell_dofs)
  (c, cache)
end

@inline function getindex!(d,a::CellDofsWithDofFixed,i::Integer)
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
    r[j] = rj
  end
  r
end
