
struct DirichletFESpace <: SingleFieldFESpace
  space::SingleFieldFESpace
end

function num_free_dofs(f::DirichletFESpace)
  num_dirichlet_dofs(f.space)
end

function zero_free_values(::Type{T},f::DirichletFESpace) where T
  zeros(T,num_free_dofs(f))
end

function get_cell_dofs(f::DirichletFESpace)
  apply(elem(-),get_cell_dofs(f.space))
end

function num_dirichlet_dofs(f::DirichletFESpace)
  num_free_dofs(f.space)
end

function zero_dirichlet_values(f::DirichletFESpace)
  T = Float64
  zero_free_values(T,f.space)
end

function num_dirichlet_tags(f::DirichletFESpace)
  1
end

function get_dirichlet_dof_tag(f::DirichletFESpace)
  ones(Int8,num_dirichlet_dofs(f))
end

function scatter_free_and_dirichlet_values(f::DirichletFESpace,fv,dv)
  scatter_free_and_dirichlet_values(f.space,dv,fv)
end

function gather_free_and_dirichlet_values(f::DirichletFESpace,cv)
  dv, fv = gather_free_and_dirichlet_values(f.space,cv)
  (fv, dv)
end

function TrialFESpace(f::DirichletFESpace)
  U = TrialFESpace(f.space)
  DirichletFESpace(U)
end

function get_cell_basis(f::DirichletFESpace)
  get_cell_basis(f.space)
end

function get_cell_dof_basis(f::DirichletFESpace)
  get_cell_dof_basis(f.space)
end

function apply_constraints_matrix_cols(f::DirichletFESpace,cm,cids)
  apply_constraints_matrix_cols(f.space,cm,cids)
end

function apply_constraints_matrix_rows(f::DirichletFESpace,cm,cids)
  apply_constraints_matrix_rows(f.space,cm,cids)
end

function apply_constraints_vector(f::DirichletFESpace,cm,cids)
  apply_constraints_vector(f.space,cm,cids)
end

