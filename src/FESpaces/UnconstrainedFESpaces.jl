
"""
Generic implementation of an unconstrained single-field FE space
Private fields and type parameters
"""
struct UnconstrainedFESpace{A,B,C} <: SingleFieldFESpace
  nfree::Int
  ndirichlet::Int
  cell_dofs::A
  cell_basis::B
  cell_dof_basis::C
  dirichlet_dof_tag::Vector{Int8}
  dirichlet_cells::Vector{Int}
  ntags::Int

  @doc """
  """
  function UnconstrainedFESpace(
    nfree::Int,
    ndirichlet::Int,
    cell_dofs::AbstractArray,
    cell_shapefuns::AbstractArray,
    cell_dof_basis::CellDofBasis,
    cell_map::AbstractArray,
    dirichlet_dof_tag::Vector{Int8},
    dirichlet_cells::Vector{Int},
    ntags) where T

    cell_basis = GenericCellBasis(cell_shapefuns,cell_map,RefStyle(cell_dof_basis))

    A = typeof(cell_dofs)
    B = typeof(cell_basis)
    C = typeof(cell_dof_basis)

    new{A,B,C}(
      nfree,
      ndirichlet,
      cell_dofs,
      cell_basis,
      cell_dof_basis,
      dirichlet_dof_tag,
      dirichlet_cells,
      ntags)
  end
end

# FESpace interface

constraint_style(::Type{<:UnconstrainedFESpace}) = Val{false}()

function num_free_dofs(f::UnconstrainedFESpace)
  f.nfree
end

function get_cell_basis(f::UnconstrainedFESpace)
  f.cell_basis
end

function zero_free_values(::Type{T},f::UnconstrainedFESpace) where T
  zeros(T,num_free_dofs(f))
end

# SingleFieldFESpace interface

function get_cell_dofs(f::UnconstrainedFESpace)
  f.cell_dofs
end

function get_cell_dof_basis(f::UnconstrainedFESpace)
  f.cell_dof_basis
end

function num_dirichlet_dofs(f::UnconstrainedFESpace)
  f.ndirichlet
end

function num_dirichlet_tags(f::UnconstrainedFESpace)
  f.ntags
end

function zero_dirichlet_values(f::UnconstrainedFESpace)
  T = Float64 # TODO
  zeros(T,num_dirichlet_dofs(f))
end

function get_dirichlet_dof_tag(f::UnconstrainedFESpace)
  f.dirichlet_dof_tag
end

function scatter_free_and_dirichlet_values(f::UnconstrainedFESpace,free_values,dirichlet_values)
  cell_dofs = get_cell_dofs(f)
  LocalToGlobalPosNegArray(cell_dofs,free_values,dirichlet_values)
end

function gather_free_and_dirichlet_values(f::UnconstrainedFESpace,cell_vals)

  cell_dofs = get_cell_dofs(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  free_vals = zero_free_values(f)
  dirichlet_vals = zero_dirichlet_values(f)
  cells = 1:length(cell_vals)

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells)

  (free_vals, dirichlet_vals)
end

function gather_dirichlet_values(f::UnconstrainedFESpace,cell_vals)

  cell_dofs = get_cell_dofs(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  free_vals = zero_free_values(f)
  dirichlet_vals = zero_dirichlet_values(f)
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
