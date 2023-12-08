
"""
Generic implementation of an unconstrained single-field FE space
Private fields and type parameters
"""
struct UnconstrainedFESpace{V,M} <: SingleFieldFESpace
  vector_type::Type{V}
  nfree::Int
  ndirichlet::Int
  cell_dofs_ids::AbstractArray
  fe_basis::CellField
  fe_dof_basis::CellDof
  cell_is_dirichlet::AbstractArray{Bool}
  dirichlet_dof_tag::Vector{Int8}
  dirichlet_cells::Vector{Int32}
  ntags::Int
  metadata::M
end

function UnconstrainedFESpace(
  vector_type::Type{V},
  nfree::Integer,
  ndirichlet::Integer,
  cell_dofs_ids::AbstractArray,
  fe_basis::CellField,
  fe_dof_basis::CellDof,
  cell_is_dirichlet::AbstractArray,
  dirichlet_dof_tag::AbstractArray,
  dirichlet_cells::AbstractArray,
  ntags::Integer) where V

  metadata = nothing
  UnconstrainedFESpace(
    V,
    nfree,
    ndirichlet,
    cell_dofs_ids,
    fe_basis,
    fe_dof_basis,
    cell_is_dirichlet,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags,
    metadata)
end

# FESpace interface

ConstraintStyle(::Type{<:UnconstrainedFESpace}) = UnConstrained()
get_free_dof_ids(f::UnconstrainedFESpace) = Base.OneTo(f.nfree)
get_fe_basis(f::UnconstrainedFESpace) = f.fe_basis
get_fe_dof_basis(f::UnconstrainedFESpace) = f.fe_dof_basis
get_cell_dof_ids(f::UnconstrainedFESpace) = f.cell_dofs_ids
get_triangulation(f::UnconstrainedFESpace) = get_triangulation(f.fe_basis)
get_dof_value_type(f::UnconstrainedFESpace{V}) where V = eltype(V)
get_vector_type(f::UnconstrainedFESpace{V}) where V = V
get_cell_is_dirichlet(f::UnconstrainedFESpace) = f.cell_is_dirichlet

# SingleFieldFESpace interface

get_dirichlet_dof_ids(f::UnconstrainedFESpace) = Base.OneTo(f.ndirichlet)
num_dirichlet_tags(f::UnconstrainedFESpace) = f.ntags
get_dirichlet_dof_tag(f::UnconstrainedFESpace) = f.dirichlet_dof_tag

function scatter_free_and_dirichlet_values(f::UnconstrainedFESpace,free_values,dirichlet_values)
  @check eltype(free_values) == eltype(dirichlet_values) """\n
  The entries stored in free_values and dirichlet_values should be of the same type.

  This error shows up e.g. when trying to build a FEFunction from a vector of integers
  if the Dirichlet values of the underlying space are of type Float64, or when the
  given free values are Float64 and the Dirichlet values ComplexF64.
  """
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(PosNegReindex(free_values,dirichlet_values)),cell_dof_ids)
end

function gather_free_and_dirichlet_values!(free_vals,dirichlet_vals,f::UnconstrainedFESpace,cell_vals)

  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  cells = 1:length(cell_vals)

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells)

  (free_vals,dirichlet_vals)
end

function gather_dirichlet_values!(dirichlet_vals,f::UnconstrainedFESpace,cell_vals)

  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  free_vals = zero_free_values(f)
  cells = f.dirichlet_cells

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells)

  dirichlet_vals
end

function  _free_and_dirichlet_values_fill!(
  free_vals,
  dirichlet_vals,
  cache_vals,
  cache_dofs,
  cell_vals,
  cell_dofs,
  cells)

  for cell in cells
    vals = getindex!(cache_vals,cell_vals,cell)
    dofs = getindex!(cache_dofs,cell_dofs,cell)
    for (i,dof) in enumerate(dofs)
      val = vals[i]
      if dof > 0
        free_vals[dof] = val
      elseif dof < 0
        dirichlet_vals[-dof] = val
      else
        @unreachable "dof ids either positive or negative, not zero"
      end
    end
  end

end
