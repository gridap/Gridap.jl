"""
    struct FESpaceWithConstantFixed{CS,RD} <: SingleFieldFESpace
      space::SingleFieldFESpace
      constraint_style::Val{CS}
      remove_dof::Val{RD}
      dof_to_remove::Int
    end
"""
struct FESpaceWithConstantFixed{CS,RD} <: SingleFieldFESpace
  space::SingleFieldFESpace
  constraint_style::Val{CS}
  remove_dof::Val{RD}
  dof_to_remove::Int
  @doc """
      FESpaceWithConstantFixed(space::SingleFieldFESpace, remove_dof::Bool,
      dof_to_remove::Int=num_free_dofs(space))
  """
  function FESpaceWithConstantFixed(space::SingleFieldFESpace,
     remove_dof::Bool,
     dof_to_remove::Int=num_free_dofs(space))
    s = "FESpaceWithConstantFixed can only be constructed from spaces without dirichlet dofs."
    @notimplementedif num_dirichlet_dofs(space) != 0 s
    @assert !remove_dof || 1 <= dof_to_remove <= num_free_dofs(space)
    cs = constraint_style(space)
    CS = get_val_parameter(cs)
    rd = Val{remove_dof}()
    RD = remove_dof
    new{CS,RD}(space,cs,rd,dof_to_remove)
  end
end

# Genuine functions
function num_free_dofs(f::FESpaceWithConstantFixed{CS,true}) where {CS}
  num_free_dofs(f.space)-1
end

function num_free_dofs(f::FESpaceWithConstantFixed{CS,false}) where {CS}
  num_free_dofs(f.space)
end

function zero_free_values(f::FESpaceWithConstantFixed)
  zeros(num_free_dofs(f))
end

function get_cell_dofs(f::FESpaceWithConstantFixed{CS,true}) where {CS}
  cell_dofs = get_cell_dofs(f.space)
  CellDofsWithDofRemoved(cell_dofs,f.dof_to_remove)
end

function get_cell_dofs(f::FESpaceWithConstantFixed{CS,false}) where {CS}
  get_cell_dofs(f.space)
end


num_dirichlet_dofs(f::FESpaceWithConstantFixed{CS,true}) where {CS} = 1

num_dirichlet_dofs(f::FESpaceWithConstantFixed{CS,false}) where {CS} = 0

function zero_dirichlet_values(f::FESpaceWithConstantFixed)
  T = Float64 # TODO
  zeros(T,num_dirichlet_dofs(f))
end


num_dirichlet_tags(f::FESpaceWithConstantFixed{CS,true}) where {CS} = 1

num_dirichlet_tags(f::FESpaceWithConstantFixed{CS,false}) where {CS} = 0

get_dirichlet_dof_tag(f::FESpaceWithConstantFixed{CS,true}) where {CS} = Int8[1,]

get_dirichlet_dof_tag(f::FESpaceWithConstantFixed{CS,false}) where {CS} = Int8[]

function scatter_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,true},fv,dv) where {CS}
  @assert length(dv) == 1
  _dv = similar(dv,eltype(dv),0)
  _fv = vcat(view(fv,1:f.dof_to_remove-1),
             dv,
             view(fv,f.dof_to_remove:length(fv))) # TODO lazy append
 scatter_free_and_dirichlet_values(f.space,_fv,_dv)
end

function scatter_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,false},fv,dv) where {CS}
  @assert length(dv) == 0
 scatter_free_and_dirichlet_values(f.space,fv,dv)
end

function gather_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,true},cv) where {CS}
  _fv, _dv = gather_free_and_dirichlet_values(f.space,cv)
  @assert length(_dv) == 0
  fv = vcat(view(_fv,1:f.dof_to_remove-1),
            view(_fv,f.dof_to_remove+1:length(_fv))) # TODO: can we avoid
                                          # allocating new memory?
  dv = view(_fv,f.dof_to_remove:f.dof_to_remove)
  (fv, dv)
end

function gather_free_and_dirichlet_values(
  f::FESpaceWithConstantFixed{CS,false},cv) where {CS}
  gather_free_and_dirichlet_values(f.space,cv)
end

function gather_free_and_dirichlet_values!(
  fv,dv,f::FESpaceWithConstantFixed{CS,true},cv) where {CS}
  _fv, _dv = gather_free_and_dirichlet_values(f.space,cv)
  @assert length(_dv) == 0
  fv .= vcat(view(_fv,1:f.dof_to_remove-1),
            view(_fv,f.dof_to_remove+1:length(_fv))) # TODO: can we avoid
                                          # allocating new memory?
  dv .= view(_fv,f.dof_to_remove:f.dof_to_remove)
  (fv, dv)
end

function gather_free_and_dirichlet_values!(
  fv,dv,f::FESpaceWithConstantFixed{CS,false},cv) where {CS}
  gather_free_and_dirichlet_values(f.space,cv)
end

function TrialFESpace(f::FESpaceWithConstantFixed{CS,RD}) where {CS,RD}
  U = TrialFESpace(f.space)
  FESpaceWithConstantFixed(U,RD,f.dof_to_remove)
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

struct CellDofsWithDofRemoved{A<:AbstractArray} <: AbstractVector{Vector{Int}}
  cell_dofs::A
  dof_to_remove::Int
end

Base.size(a::CellDofsWithDofRemoved) = (length(a.cell_dofs),)

Base.IndexStyle(::Type{<:CellDofsWithDofRemoved}) = IndexLinear()

function Base.getindex(a::CellDofsWithDofRemoved,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function array_cache(a::CellDofsWithDofRemoved)
  @assert eltype(a.cell_dofs) == Vector{Int}
  b = testitem(a.cell_dofs)
  c = CachedArray(b)
  cache = array_cache(a.cell_dofs)
  (c, cache)
end

@inline function getindex!(d,a::CellDofsWithDofRemoved,i::Integer)
  c, cache = d
  b = getindex!(cache,a.cell_dofs,i)
  setsize!(c,size(b))
  r = c.array
  for j in 1:length(b)
    bj = b[j]
    if bj == a.dof_to_remove
      rj = -1
    elseif bj < a.dof_to_remove
      rj = bj
    else
      rj = bj-1
    end
    r[j] = rj
  end
  r
end
