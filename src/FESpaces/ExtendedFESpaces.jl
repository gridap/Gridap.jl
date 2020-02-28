
# TODO This is temporary
using Gridap.Arrays: Reindexed

"""
"""
struct ExtendedFESpace <: SingleFieldFESpace
  space::SingleFieldFESpace
  trian::RestrictedTriangulation
end

constraint_style(::Type{ExtendedFESpace}) = Val{false}()

function get_cell_dofs(f::ExtendedFESpace)
  Reindexed(get_cell_dofs(f.space),f.trian.oldcell_to_cell)
end

function get_cell_basis(f::ExtendedFESpace)
  cell_basis = get_cell_basis(f.space)
  a = get_array(cell_basis)
  array = Reindexed(a,f.trian.oldcell_to_cell)
  b = get_cell_map(cell_basis)
  cm = get_cell_map(f.trian.oldtrian)
  trial_style = TrialStyle(cell_basis)
  GenericCellBasis(trial_style,array,cm)
end

function get_cell_dof_basis(f::ExtendedFESpace)
  Reindexed(get_cell_dof_basis(f.space),f.trian.oldcell_to_cell)
end

function scatter_free_and_dirichlet_values(f::ExtendedFESpace,fv,dv)
  a = scatter_free_and_dirichlet_values(f.space,fv,dv)
  Reindexed(a,f.trian.oldcell_to_cell)
end

function gather_free_and_dirichlet_values(f::ExtendedFESpace,cv)
  _cv = reindex(cv,f.trian.cell_to_oldcell)
  gather_free_and_dirichlet_values(f.space,_cv)
end

function num_free_dofs(f::ExtendedFESpace)
  num_free_dofs(f.space)
end

function zero_free_values(::Type{T},f::ExtendedFESpace) where T
  zeros(T,num_free_dofs(f))
end

function num_dirichlet_dofs(f::ExtendedFESpace)
  num_dirichlet_dofs(f.space)
end

function zero_dirichlet_values(f::ExtendedFESpace)
  zero_dirichlet_values(f.space)
end

function num_dirichlet_tags(f::ExtendedFESpace)
  num_dirichlet_tags(f.space)
end

function get_dirichlet_dof_tag(f::ExtendedFESpace)
  get_dirichlet_dof_tag(f.space)
end

function TrialFESpace(f::ExtendedFESpace)
  U = TrialFESpace(f.space)
  ExtendedFESpace(U,f.trian)
end

