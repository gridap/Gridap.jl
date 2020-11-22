
"""
Generic implementation of an unconstrained single-field FE space
Private fields and type parameters
"""
struct UnconstrainedFESpace{V} <: SingleFieldFESpace
  vector_type::Type{V}
  nfree::Int
  ndirichlet::Int
  cell_dofs_ids::AbstractArray
  cell_shapefuns::CellField
  cell_dof_basis::CellDof
  dirichlet_dof_tag::Vector{Int8}
  dirichlet_cells::Vector{Int32}
  ntags::Int
end

# FESpace interface

ConstraintStyle(::Type{<:UnconstrainedFESpace}) = UnConstrained()
num_free_dofs(f::UnconstrainedFESpace) = f.nfree
zero_free_values(f::UnconstrainedFESpace) = allocate_vector(f.vector_type,num_free_dofs(f))
get_cell_shapefuns(f::UnconstrainedFESpace) = f.cell_shapefuns
get_cell_dof_basis(f::UnconstrainedFESpace) = f.cell_dof_basis
get_cell_dof_ids(f::UnconstrainedFESpace) = f.cell_dofs_ids
get_triangulation(f::UnconstrainedFESpace) = get_triangulation(f.cell_shapefuns)
get_dof_value_type(f::UnconstrainedFESpace{V}) where V = eltype(V)
get_vector_type(f::UnconstrainedFESpace{V}) where V = V

# SingleFieldFESpace interface

num_dirichlet_dofs(f::UnconstrainedFESpace) = f.ndirichlet
num_dirichlet_tags(f::UnconstrainedFESpace) = f.ntags
zero_dirichlet_values(f::UnconstrainedFESpace) = allocate_vector(f.vector_type,num_dirichlet_dofs(f))
get_dirichlet_dof_tag(f::UnconstrainedFESpace) = f.dirichlet_dof_tag

function scatter_free_and_dirichlet_values(f::UnconstrainedFESpace,free_values,dirichlet_values)
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
        @unreachable "dof ids either positibe or negative, not zero"
      end
    end
  end

end
