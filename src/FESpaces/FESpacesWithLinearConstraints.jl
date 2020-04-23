
# Notation in this file
# dof: a dof (either free or Dirichlet) in the original space
# fdof: a free dof in the original space
# ddof: a Dirichlet dof in the original space
# mdof: a master dof (i.e. a free dof in the constrained space)
# mddof: Either a mdof or a ddof (the sign of the value tells if it is one or the other)
# ldof: a local dof in a cell
# ledof: local "extended" dof in a cell
#


struct FESpaceWithLinearConstraints <: SingleFieldFESpace
  space::SingleFieldFESpace
  num_mdofs::Int
  mdof_to_fdof::Vector{Int}
  fdof_to_mddofs::Table{Int,Int32}
  dof_to_coeffs::Table{Int,Int32}
  cell_to_ledof_to_mddof::Table{Int,Int32}
end

function FESpaceWithLinearConstraints(
  space::SingleFieldFESpace, dof_to_dofs::Table, dof_to_coeffs::Table)

  # Create master dof (mdof) numeration
  fdof_to_mdof, num_mdofs = _setup_fdof_to_mdof(dof_to_dofs)

  cell_to_ledof_to_mddof = _setup_cell_to_ldof_to_mddof(dof_to_dofs,fdof_to_mdof)


end

function num_free_dofs(f::FESpaceWithLinearConstraints)
  f.num_mdofs
end

function zero_free_values(::Type{T},f::FESpaceWithLinearConstraints) where T
  zeros(T,num_free_dofs(f))
end

function get_cell_dofs(f::FESpaceWithLinearConstraints)
  f.cell_to_ledof_to_mddof
end

num_dirichlet_dofs(f::FESpaceWithLinearConstraints) = num_dirichlet_dofs(f.space)

num_dirichlet_tags(f::FESpaceWithLinearConstraints) = num_dirichlet_tags(f.space)

get_dirichlet_dof_tag(f::FESpaceWithLinearConstraints) = get_dirichlet_dof_tag(f.space)

function scatter_free_and_dirichlet_values(f::FESpaceWithLinearConstraints,mvals,dvals)
  mdof_to_val = mvals
  ddof_to_val = dvals
  fdof_to_mddofs = f.fdof_to_mddofs
  fdof_to_coeffs = f.dof_to_coeffs
  fdof_to_val = zero_free_values(eltype(mvals),f.space)
  _setup_fdof_to_val!(fdof_to_val,mdof_to_val,ddof_to_val,fdof_to_mddofs,fdof_to_coeffs)
  scatter_free_and_dirichlet_values(f.space,fdof_to_val,ddof_to_val)
end

function gather_free_and_dirichlet_values(f::FESpaceWithLinearConstraints,cv)
  fdof_to_val, ddof_to_val = gather_free_and_dirichlet_values(f.space,cv)
  mdof_to_val = fdof_to_val[f.mdof_to_fdof]
  mdof_to_val, ddof_to_val
end

function get_cell_basis(f::FESpaceWithLinearConstraints)
  get_cell_basis(f.space)
end

function get_cell_dof_basis(f::FESpaceWithLinearConstraints)
  get_cell_dof_basis(f.space)
end

constraint_style(::Type{<:FESpaceWithLinearConstraints}) = Val{true}()

# Helpers

function _setup_fdof_to_val!(fdof_to_val,mdof_to_val,ddof_to_val,fdof_to_mddofs,fdof_to_coeffs)
  for fdof in 1:length(fdof_to_val)
    pini = fdof_to_mddofs.ptrs[fdof]
    pend = fdof_to_mddofs.ptrs[fdof+1]-1
    for p in pini:pend
      mddof = fdof_to_mddofs.data[p]
      coeff = fdof_to_coeffs.data[p]
      if mddof > 0
        mdof = mddof
        val = mdof_to_val[mdof]
      else
        ddof = -mddof
        val = ddof_to_val[ddof]
      end
      fdof_to_val[fdof] += val*coeff
    end
  end
end


