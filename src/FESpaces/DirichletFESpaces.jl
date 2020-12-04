
"""
    struct DirichletFESpace <: SingleFieldFESpace
      space::SingleFieldFESpace
    end
"""
struct DirichletFESpace{S<:SingleFieldFESpace} <: SingleFieldFESpace
  space::S
  function DirichletFESpace(space::SingleFieldFESpace)
    new{typeof(space)}(space)
  end
end

ConstraintStyle(::Type{DirichletFESpace{B}}) where B = ConstraintStyle(B)

CellField(t::DirichletFESpace,cell_vals) = CellField(t.space,cell_vals)

get_cell_isconstrained(f::DirichletFESpace) = get_cell_isconstrained(f.space)

get_cell_constraints(f::DirichletFESpace) = get_cell_constraints(f.space)

function num_free_dofs(f::DirichletFESpace)
  num_dirichlet_dofs(f.space)
end

function get_vector_type(f::DirichletFESpace)
  get_vector_type(f.space)
end

function get_cell_dof_ids(f::DirichletFESpace)
  lazy_map(Broadcasting(-),get_cell_dof_ids(f.space))
end

function num_dirichlet_dofs(f::DirichletFESpace)
  num_free_dofs(f.space)
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

function get_cell_shapefuns(f::DirichletFESpace)
  get_cell_shapefuns(f.space)
end

function get_cell_shapefuns_trial(f::DirichletFESpace)
  get_cell_shapefuns_trial(f.space)
end

function get_cell_dof_basis(f::DirichletFESpace)
  get_cell_dof_basis(f.space)
end

function get_triangulation(f::DirichletFESpace)
  get_triangulation(f.space)
end
