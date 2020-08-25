
struct TrialFESpace{B} <: SingleFieldFESpace
  space::SingleFieldFESpace
  dirichlet_values::AbstractVector
  cell_basis::CellField
  constraint_style::Val{B}

  function TrialFESpace(dirichlet_values::AbstractVector,space::SingleFieldFESpace)
    cell_basis = trialize_cell_basis(get_cell_basis(space))
    cs = constraint_style(space)
    B = get_val_parameter(cs)
    new{B}(space,dirichlet_values,cell_basis,cs)
  end
end

"""
"""
function TrialFESpace(space::SingleFieldFESpace)
  dirichlet_values = get_dirichlet_values(space)
  TrialFESpace(dirichlet_values,space)
end

"""
"""
function TrialFESpace(space::SingleFieldFESpace,objects)
  dirichlet_values = compute_dirichlet_values_for_tags(space,objects)
  TrialFESpace(dirichlet_values,space)
end

"""
"""
function TrialFESpace!(dir_values::AbstractVector,space::SingleFieldFESpace,objects)
  dir_values = compute_dirichlet_values_for_tags!(dir_values,space,objects)
  TrialFESpace(dir_values,space)
end

"""
"""
function TrialFESpace!(space::TrialFESpace,objects)
  dir_values = get_dirichlet_values(space)
  dir_values = compute_dirichlet_values_for_tags!(dir_values,space,objects)
  space
end

function TrialFESpace(space::TrialFESpace)
  space
end

function HomogeneousTrialFESpace(U::FESpace)
  dirichlet_values = zero_dirichlet_values(U)
  TrialFESpace(dirichlet_values,U)
end

function HomogeneousTrialFESpace!(dirichlet_values::AbstractVector,U::FESpace)
  fill!(dirichlet_values,zero(eltype(dirichlet_values)))
  TrialFESpace(dirichlet_values,U)
end

# Genuine functions

get_dirichlet_values(f::TrialFESpace) = f.dirichlet_values

get_cell_basis(f::TrialFESpace) = f.cell_basis

# Delegated functions

get_cell_axes(t::TrialFESpace)= get_cell_axes(t.space)

get_cell_axes_with_constraints(t::TrialFESpace)= get_cell_axes_with_constraints(t.space)

CellData.CellField(t::TrialFESpace,cell_vals) = CellField(t.space,cell_vals)

constraint_style(::Type{<:TrialFESpace{B}}) where B = Val{B}()

get_cell_dof_basis(f::TrialFESpace) = get_cell_dof_basis(f.space)

num_free_dofs(f::TrialFESpace) = num_free_dofs(f.space)

zero_free_values(f::TrialFESpace) = zero_free_values(f.space)

get_cell_dofs(f::TrialFESpace) = get_cell_dofs(f.space)

num_dirichlet_dofs(f::TrialFESpace) = num_dirichlet_dofs(f.space)

zero_dirichlet_values(f::TrialFESpace) = zero_dirichlet_values(f.space)

num_dirichlet_tags(f::TrialFESpace) = num_dirichlet_tags(f.space)

get_dirichlet_dof_tag(f::TrialFESpace) = get_dirichlet_dof_tag(f.space)

scatter_free_and_dirichlet_values(f::TrialFESpace,fv,dv) = scatter_free_and_dirichlet_values(f.space,fv,dv)

gather_free_and_dirichlet_values(f::TrialFESpace,cv) = gather_free_and_dirichlet_values(f.space,cv)

gather_free_and_dirichlet_values!(fv,dv,f::TrialFESpace,cv) = gather_free_and_dirichlet_values!(fv,dv,f.space,cv)

gather_dirichlet_values(f::TrialFESpace,cv) = gather_dirichlet_values(f.space,cv)

gather_dirichlet_values!(dv,f::TrialFESpace,cv) = gather_dirichlet_values!(dv,f.space,cv)

gather_free_values(f::TrialFESpace,cv) = gather_free_values(f.space,cv)

gather_free_values!(fv,f::TrialFESpace,cv) = gather_free_values!(fv,f.space,cv)

get_cell_isconstrained(f::TrialFESpace) = get_cell_isconstrained(f.space)

get_cell_constraints(f::TrialFESpace) = get_cell_constraints(f.space)

