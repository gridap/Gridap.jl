
"""
    struct FESpaceWithLastDofRemoved <: SingleFieldFESpace
      space::SingleFieldFESpace
    end
"""
struct FESpaceWithLastDofRemoved{B} <: SingleFieldFESpace
  space::SingleFieldFESpace
  constraint_style::Val{B}
  @doc """
      FESpaceWithLastDofRemoved(space::SingleFieldFESpace)
  """
  function FESpaceWithLastDofRemoved(space::SingleFieldFESpace)
    s = "FESpaceWithLastDofRemoved can only be constructed from spaces without dirichlet dofs."
    @notimplementedif num_dirichlet_dofs(space) != 0 s
    cs = constraint_style(space)
    B = get_val_parameter(cs)
    new{B}(space,cs)
  end
end

# Genuine functions

function num_free_dofs(f::FESpaceWithLastDofRemoved)
  num_free_dofs(f.space) - 1
end

function zero_free_values(::Type{T},f::FESpaceWithLastDofRemoved) where T
  zeros(T,num_free_dofs(f))
end

function get_cell_dofs(f::FESpaceWithLastDofRemoved)
  cell_dofs = get_cell_dofs(f.space)
  n = num_free_dofs(f)
  CellDofsWithLastRemoved(cell_dofs,n)
end

num_dirichlet_dofs(f::FESpaceWithLastDofRemoved) = 1

function zero_dirichlet_values(f::FESpaceWithLastDofRemoved)
  T = Float64 # TODO
  zeros(T,num_dirichlet_dofs(f))
end

num_dirichlet_tags(f::FESpaceWithLastDofRemoved) = 1

get_dirichlet_dof_tag(f::FESpaceWithLastDofRemoved) = Int8[1,]

function scatter_free_and_dirichlet_values(
  f::FESpaceWithLastDofRemoved,fv,dv)
  @assert length(dv) == 1
  _dv = similar(dv,eltype(dv),0)
  _fv = vcat(fv,dv) # TODO lazy append
 scatter_free_and_dirichlet_values(f.space,_fv,_dv)
end

function gather_free_and_dirichlet_values(
  f::FESpaceWithLastDofRemoved,cv)
  _fv, _dv = gather_free_and_dirichlet_values(f.space,cv)
  @assert length(_dv) == 0
  l = length(_fv)
  fv = SubVector(_fv,1,l-1)
  dv = SubVector(_fv,l,l)
  (fv, dv)
end

function TrialFESpace(f::FESpaceWithLastDofRemoved)
  U = TrialFESpace(f.space)
  FESpaceWithLastDofRemoved(U)
end

# Delegated functions

function get_cell_basis(f::FESpaceWithLastDofRemoved)
  get_cell_basis(f.space)
end

function get_cell_dof_basis(f::FESpaceWithLastDofRemoved)
  get_cell_dof_basis(f.space)
end

constraint_style(::Type{FESpaceWithLastDofRemoved{B}}) where B = Val{B}()

function get_constraint_kernel_matrix_cols(f::FESpaceWithLastDofRemoved)
  get_constraint_kernel_matrix_cols(f.space)
end

function get_constraint_kernel_matrix_rows(f::FESpaceWithLastDofRemoved)
  get_constraint_kernel_matrix_rows(f.space)
end

function get_constraint_kernel_vector(f::FESpaceWithLastDofRemoved)
  get_constraint_kernel_vector(f.space)
end

# Helpers

struct CellDofsWithLastRemoved{A<:AbstractArray} <: AbstractVector{Vector{Int}}
  cell_dofs::A
  n::Int
end

Base.size(a::CellDofsWithLastRemoved) = (length(a.cell_dofs),)

Base.IndexStyle(::Type{<:CellDofsWithLastRemoved}) = IndexLinear()

function Base.getindex(a::CellDofsWithLastRemoved,i::Integer)
  cache = array_cache(a)
  getindex!(cache,a,i)
end

function array_cache(a::CellDofsWithLastRemoved)
  @assert eltype(a.cell_dofs) == Vector{Int}
  b = testitem(a.cell_dofs)
  c = CachedArray(b)
  cache = array_cache(a.cell_dofs)
  (c, cache)
end

@inline function getindex!(d,a::CellDofsWithLastRemoved,i::Integer)
  c, cache = d
  b = getindex!(cache,a.cell_dofs,i)
  setsize!(c,size(b))
  r = c.array
  for j in 1:length(b)
    bj = b[j]
    if bj > a.n
      rj = -1
    else
      rj = bj
    end
    r[j] = rj
  end
  r
end

