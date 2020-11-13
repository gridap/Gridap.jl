
struct TrialFESpace{S} <: SingleFieldFESpace
  dirichlet_values::AbstractVector
  space::S
  function TrialFESpace(dirichlet_values::AbstractVector,space::SingleFieldFESpace)
    new{typeof(space)}(dirichlet_values,space)
  end
end

"""
"""
function TrialFESpace(space::SingleFieldFESpace)
  space
end

#function TrialFESpace(space::TrialFESpace)
#  space
#end

"""
"""
function TrialFESpace(space::SingleFieldFESpace,objects)
  dirichlet_values = compute_dirichlet_values_for_tags(space,objects)
  TrialFESpace(dirichlet_values,space)
end

"""
"""
function TrialFESpace!(dir_values::AbstractVector,space::SingleFieldFESpace,objects)
  dir_values_scratch = zero_dirichlet_values(space)
  dir_values = compute_dirichlet_values_for_tags!(dir_values,dir_values_scratch,space,objects)
  TrialFESpace(dir_values,space)
end

"""
"""
function TrialFESpace!(space::TrialFESpace,objects)
  dir_values = get_dirichlet_values(space)
  dir_values_scratch = zero_dirichlet_values(space)
  dir_values = compute_dirichlet_values_for_tags!(dir_values,dir_values_scratch,space,objects)
  space
end

# Allow do-block syntax

function TrialFESpace(f::Function,space::SingleFieldFESpace)
  TrialFESpace(space,f)
end

function TrialFESpace!(f::Function,dir_values::AbstractVector,space::SingleFieldFESpace)
  TrialFESpace!(dir_values,space,f)
end

function TrialFESpace!(f::Function,space::TrialFESpace)
  TrialFESpace!(space,f)
end

# Remove Dirichlet from the given space

function HomogeneousTrialFESpace(U::SingleFieldFESpace)
  dirichlet_values = zero_dirichlet_values(U)
  TrialFESpace(dirichlet_values,U)
end

function HomogeneousTrialFESpace!(dirichlet_values::AbstractVector,U::SingleFieldFESpace)
  fill!(dirichlet_values,zero(eltype(dirichlet_values)))
  TrialFESpace(dirichlet_values,U)
end

# Genuine functions

get_dirichlet_values(f::TrialFESpace) = f.dirichlet_values

# Delegated functions

num_free_dofs(f::TrialFESpace) = num_free_dofs(f.space)

zero_free_values(f::TrialFESpace) = zero_free_values(f.space)

get_triangulation(f::TrialFESpace) = get_triangulation(f.space)

get_dof_value_type(f::TrialFESpace) = get_dof_value_type(f.space)

get_vector_type(f::TrialFESpace) = get_vector_type(f.space)

get_cell_dof_ids(f::TrialFESpace) = get_cell_dof_ids(f.space)

get_cell_shapefuns(f::TrialFESpace) = get_cell_shapefuns(f.space)

get_cell_shapefuns_trial(f::TrialFESpace) = get_cell_shapefuns_trial(f.space)

get_cell_dof_basis(f::TrialFESpace) = get_cell_dof_basis(f.space)

ConstraintStyle(::Type{<:TrialFESpace{B}}) where B = ConstraintStyle(B)

get_cell_isconstrained(f::TrialFESpace) = get_cell_isconstrained(f.space)

get_cell_constraints(f::TrialFESpace) = get_cell_constraints(f.space)

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

