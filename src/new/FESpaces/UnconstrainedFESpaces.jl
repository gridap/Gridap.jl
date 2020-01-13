
"""
Generic implementation of an unconstrained single-field FE space
Private fields and type parameters
"""
struct UnsconstrainedFESpace{T,A,B,C} <: SingleFieldFESpace
  nfree::Int
  ndirichlet::Int
  cell_dofs::A
  cell_fe_basis::B
  cell_dof_basis::C
  dirichlet_dof_tag::Vector{Int8}

  @doc """
  """
  function UnsconstrainedFESpace(
    ::Type{T},
    nfree::Int,
    ndirichlet::Int,
    cell_dofs::AbstractArray,
    cell_shapefuns::AbstractArray,
    cell_dof_basis::AbstractArray,
    cell_map::AbstractArray,
    dirichlet_dof_tag::Vector{Int8}) where T

    cell_fe_basis = CellShapeFunsWithMap(cell_shapefuns,cell_map)

    A = typeof(cell_dofs)
    B = typeof(cell_fe_basis)
    C = typeof(cell_dof_basis)

    new{T,A,B,C}(
      nfree,
      ndirichlet,
      cell_dofs,
      cell_fe_basis,
      cell_dof_basis,
      dirichlet_dof_tag)
  end
end

# FESpace interface

function num_free_dofs(f::UnsconstrainedFESpace)
  f.nfree
end

function get_cell_fe_basis(f::UnsconstrainedFESpace)
  f.cell_fe_basis
end

function zero_free_values(f::UnsconstrainedFESpace{T}) where T
  zeros(T,num_free_dofs(f))
end

function apply_constraints_matrix_cols(f::UnsconstrainedFESpace,cellmat,cellids)
  cellmat
end

function apply_constraints_matrix_rows(f::UnsconstrainedFESpace,cellmat,cellids)
  cellmat
end

function apply_constraints_vector(f::UnsconstrainedFESpace,cellvec,cellids)
  cellvec
end

# SingleFieldFESpace interface

function get_cell_dofs(f::UnsconstrainedFESpace)
  f.cell_dofs
end

function get_cell_dof_basis(f::UnsconstrainedFESpace)
  f.cell_dof_basis
end

function num_dirichlet_dofs(f::UnsconstrainedFESpace)
  f.ndirichlet
end

function zero_dirichlet_values(f::UnsconstrainedFESpace{T}) where T
  zeros(T,num_dirichlet_dofs(f))
end

function get_dirichlet_dof_tag(f::UnsconstrainedFESpace)
  f.dirichlet_dof_tag
end

function scatter_free_and_dirichlet_values(f::UnsconstrainedFESpace,free_values,dirichlet_values)
  cell_dofs = get_cell_dofs(f)
  LocalToGlobalPosNegArray(cell_dofs,free_vals,dirichlet_values)
end

function gather_free_and_dirichlet_values(f::UnsconstrainedFESpace,cell_vals)

  cell_dofs = get_cell_dofs(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  free_vals = zero_free_values(f)
  dirichlet_vals = zero_dirichlet_values(f)

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs)

  (free_vals, dirichlet_vals)
end

function  _free_and_dirichlet_values_fill!(
  free_vals,
  dirichlet_vals,
  cache_vals,
  cache_dofs,
  cell_vals,
  cell_dofs)

  for cell in 1:length(cell_vals)
    vals = getindex!(cache_vals,cell_vals)
    dofs = getindex!(cache_dofs,cell_dofs)
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
