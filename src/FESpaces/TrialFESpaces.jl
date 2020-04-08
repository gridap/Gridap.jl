
struct TrialFESpace{B} <: SingleFieldFESpace
  space::SingleFieldFESpace
  dirichlet_values::AbstractVector
  cell_basis::CellBasis
  constraint_style::Val{B}
end

"""
"""
function TrialFESpace(space::SingleFieldFESpace)
  dirichlet_values = get_dirichlet_values(space)
  cell_basis = _prepare_trial_cell_basis(space)
  TrialFESpace(space,dirichlet_values,cell_basis,constraint_style(space))
end

"""
"""
function TrialFESpace(space::SingleFieldFESpace,objects)
  dirichlet_values = compute_dirichlet_values_for_tags(space,objects)
  cell_basis = _prepare_trial_cell_basis(space)
  TrialFESpace(space,dirichlet_values,cell_basis,constraint_style(space))
end

"""
"""
function TrialFESpace(dir_values::AbstractVector,space::SingleFieldFESpace,objects)
  dir_values = compute_dirichlet_values_for_tags!(dir_values,space,objects)
  cell_basis = _prepare_trial_cell_basis(space)
  TrialFESpace(space,dir_values,cell_basis,constraint_style(space))
end

"""
"""
function TrialFESpace!(space::SingleFieldFESpace,objects)
  dir_values = get_dirichlet_values(space)
  dir_values = compute_dirichlet_values_for_tags!(dir_values,space,objects)
  space
end

function TrialFESpace(space::TrialFESpace)
  space
end

function  _prepare_trial_cell_basis(space)
  cb = get_cell_basis(space)
  a = get_array(cb)
  cm = get_cell_map(cb)
  trial_style = Val{true}()
  cell_basis = GenericCellBasis(trial_style,a,cm,RefStyle(cb))
end

# Genuine functions

get_dirichlet_values(f::TrialFESpace) = f.dirichlet_values

get_cell_basis(f::TrialFESpace) = f.cell_basis

# Delegated functions

constraint_style(::Type{<:TrialFESpace{B}}) where B = Val{B}()

get_cell_dof_basis(f::TrialFESpace) = get_cell_dof_basis(f.space)

num_free_dofs(f::TrialFESpace) = num_free_dofs(f.space)

zero_free_values(::Type{T},f::TrialFESpace) where T = zero_free_values(T,f.space)

get_cell_dofs(f::TrialFESpace) = get_cell_dofs(f.space)

num_dirichlet_dofs(f::TrialFESpace) = num_dirichlet_dofs(f.space)

zero_dirichlet_values(f::TrialFESpace) = zero_dirichlet_values(f.space)

num_dirichlet_tags(f::TrialFESpace) = num_dirichlet_tags(f.space)

get_dirichlet_dof_tag(f::TrialFESpace) = get_dirichlet_dof_tag(f.space)

scatter_free_and_dirichlet_values(f::TrialFESpace,fv,dv) = scatter_free_and_dirichlet_values(f.space,fv,dv)

gather_free_and_dirichlet_values(f::TrialFESpace,cv) = gather_free_and_dirichlet_values(f.space,cv)

gather_dirichlet_values(f::TrialFESpace,cv) = gather_dirichlet_values(f.space,cv)

gather_free_values(f::TrialFESpace,cv) = gather_free_values(f.space,cv)
