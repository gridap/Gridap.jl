
"""
    struct DirichletFESpace <: SingleFieldFESpace
      space::SingleFieldFESpace
    end
"""
struct DirichletFESpace{B} <: SingleFieldFESpace
  space::SingleFieldFESpace
  constraint_style::Val{B}
  function DirichletFESpace(space::SingleFieldFESpace)
    cs = constraint_style(space)
    B = get_val_parameter(cs)
    new{B}(space,cs)
  end
end

constraint_style(::Type{DirichletFESpace{B}}) where B = Val{B}()

function get_constraint_kernel_matrix_cols(f::DirichletFESpace)
  get_constraint_kernel_matrix_cols(f.space)
end

function get_constraint_kernel_matrix_rows(f::DirichletFESpace)
  get_constraint_kernel_matrix_rows(f.space)
end

function get_constraint_kernel_vector(f::DirichletFESpace)
  get_constraint_kernel_vector(f.space)
end

function num_free_dofs(f::DirichletFESpace)
  num_dirichlet_dofs(f.space)
end

function zero_free_values(f::DirichletFESpace)
  zero_dirichlet_values(f.space)
end

function get_cell_dofs(f::DirichletFESpace)
  apply(elem(-),get_cell_dofs(f.space))
end

function num_dirichlet_dofs(f::DirichletFESpace)
  num_free_dofs(f.space)
end

function zero_dirichlet_values(f::DirichletFESpace)
  zero_free_values(f.space)
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

function gather_free_and_dirichlet_values!(fv,dv,f::DirichletFESpace,cv)
  gather_free_and_dirichlet_values!(dv,fv,f.space,cv)
  (fv,dv)
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

