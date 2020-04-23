
# Background:
#
# We build a novel fe space form a given fe space plus a set of linear constraints.
# We accept any single field fe space as input.  In particular, the given fe space can also be defied
# via constraints.
#
# Assumptions:
#
#  - A constrained dof depends only on master dofs, i.e., only one level of constraints
#    This can be easily relaxed in the future.
#
#  - We can add constraints to Dirichlet dofs but their masters have to be Dirichlet as well.
#    However, free dofs can be constrained both by free and Dirichlet.
# 
# Notation in this file:
#
#  - dof: a dof (either free or Dirichlet) in the original space.
#    The sign of the value tells if it is one or the other.
#
#  - fdof: a free dof in the original space. Always positive.
#
#  - ddof: a Dirichlet dof in the original space.
#    Always positive, that is ddof = -dof for its corresponding dof.
#
#  - mdof: a master dof (i.e. a dof in the constrained space).
#    It can be free or Dirichlet. The sign of the value tells if it is one or the other.
#
#  - fmdof: a free master dof. Always positive.
#
#  - dmdof: a Dirichlet master dof. Always positive, that is dmdof = -mdof for its corresponding mdof.
#
#  - ldof: a local dof in a cell for the original space. Always positive.
#
#  - lmdof: local master dof in a cell. That is, a local dofs in a cell for the constrained space.
#    Always positive.
#
#  - ludof: A local dof in the unconstrained interpolation space of the underlying cell. Always positive.
#
#
#  Note that only values that are marked as "Always positive" can be used to index arrays.


struct FESpaceWithLinearConstraints <: SingleFieldFESpace
  space::SingleFieldFESpace
  fmdof_to_fdof::Vector
  dmdof_to_ddof::Vector
  fdof_to_mdofs::Table
  fdof_to_coeffs::Table
  ddof_to_dmdofs::Table
  ddof_to_coeffs::Table
  cell_to_lmdof_to_mdof::Table
  cell_to_lmdof_to_ldofs_and_coeffs

  function FESpaceWithLinearConstraints(
    space::SingleFieldFESpace,
    fdof_to_dofs::Table,
    fdof_to_coeffs::Table,
    ddof_to_ddofs::Table,
    ddof_to_coeffs::Table)
  
    fmdof_to_fdof, dmdof_to_ddof = _find_master_dofs(fdof_to_dofs,ddof_to_ddofs)

    fdof_to_mdofs, ddof_to_dmdofs = _renumber_constraints(
      fdof_to_dofs,ddof_to_ddofs,fmdof_to_fdof, dmdof_to_ddof)
  
    cell_to_ldof_to_dof = Table(get_cell_dofs(space))

    cell_to_lmdof_to_mdof = _setup_cell_to_lmdof_to_mdof(
      cell_to_ldof_to_dof,fdof_to_mdofs,ddof_to_dmdofs)
  
  
  end

end

function _setup_cell_to_lmdof_to_mdof(cell_to_ldof_to_dof,fdof_to_mdofs,ddof_to_dmdofs)

  n_cells = length(cell_to_ldof_to_dof)
  cell_to_lmdof_to_mdof_ptrs = zeros(eltype(fdof_to_mdofs.ptrs),n_cells)

  for cell in 1:n_cells
    mdofs = Set{Int}()
    pini = cell_to_ldof_to_dof.ptrs[cell]
    pend = cell_to_ldof_to_dof.ptrs[cell+1]-1
    for p in pini:pend
      dof = cell_to_ldof_to_dof.data[p]
      if dof>0
        fdof = dof
        qini = fdof_to_mdofs.ptrs[fdof]
        qend = fdof_to_mdofs.ptrs[fdof+1]-1
        for q in qini:qend
          mdof = fdof_to_mdofs.data[q]
          push!(mdofs,mdof)
        end
      else
        ddof = -dof
        qini = ddof_to_dmdofs.ptrs[ddof]
        qend = ddof_to_dmdofs.ptrs[ddof+1]-1
        for q in qini:qend
          dmdof = ddof_to_dmdofs.data[q]
          mdof = -dmdof
          push!(mdofs,mdof)
        end
      end
    end
    cell_to_lmdof_to_mdof_ptrs[cell+1] = length(mdofs)
  end

  length_to_ptrs!(cell_to_lmdof_to_mdof_ptrs)
  ndata = cell_to_lmdof_to_mdof_ptrs[end]-1

  cell_to_lmdof_to_mdof_data = zeros(eltype(fdof_to_mdofs.data),ndata)

  for cell in 1:n_cells
    mdofs = Set{Int}()
    pini = cell_to_ldof_to_dof.ptrs[cell]
    pend = cell_to_ldof_to_dof.ptrs[cell+1]-1
    for p in pini:pend
      dof = cell_to_ldof_to_dof.data[p]
      if dof>0
        fdof = dof
        qini = fdof_to_mdofs.ptrs[fdof]
        qend = fdof_to_mdofs.ptrs[fdof+1]-1
        for q in qini:qend
          mdof = fdof_to_mdofs.data[q]
          push!(mdofs,mdof)
        end
      else
        ddof = -dof
        qini = ddof_to_dmdofs.ptrs[ddof]
        qend = ddof_to_dmdofs.ptrs[ddof+1]-1
        for q in qini:qend
          dmdof = ddof_to_dmdofs.data[q]
          mdof = -dmdof
          push!(mdofs,mdof)
        end
      end
    end

    o = cell_to_lmdof_to_mdof_ptrs[cell]-1
    for (lmdof, mdof) in enumerate(mdofs)
      cell_to_lmdof_to_mdof_data[o+lmdof] = mdof
    end
  end

  Table(cell_to_lmdof_to_mdof_data,cell_to_lmdof_to_mdof_ptrs)
end

function _renumber_constraints(fdof_to_dofs,ddof_to_ddofs,fmdof_to_fdof,dmdof_to_ddof)

  fdof_to_fmdof = zeros(eltype(fmdof_to_fdof),length(fdof_to_dofs))
  fdof_to_fmdof .= 1:length(fmdof_to_fdof)

  ddof_to_dmdof = zeros(eltype(dmdof_to_ddof),length(ddof_to_ddofs))
  ddof_to_dmdof .= 1:length(dmdof_to_ddof)

  fdof_to_mdofs = Table(LocalToGlobalPosNegArray(fdof_to_dofs,fdof_to_fmdof,apply(-,ddof_to_dmdof)))
  ddof_to_dmdofs = Table(LocalToGlobalArray(ddof_to_ddofs,ddof_to_dmdof))

  fdof_to_mdofs, ddof_to_dmdofs
end

function _find_master_dofs(fdof_to_dofs,ddof_to_ddofs)

  n_fdofs = length(fdof_to_dofs)
  n_ddofs = length(ddof_to_ddofs)

  fdof_to_ismaster = fill(false,n_fdofs)
  ddof_to_ismaster = fill(false,n_ddofs)

  for fdof in 1:n_fdofs
    pini = fdof_to_dofs.ptrs[fdof]
    pend = fdof_to_dofs.ptrs[fdof+1]-1
    for p in pini:pend
      dof = fdof_to_dofs.data[p]
      if dof > 0
        @assert (fdof_to_dofs.ptrs[dof+1]-fdof_to_dofs.ptrs[dof]) == 1
        @assert fdof_to_dofs.data[fdof_to_dofs.ptrs[dof]] == dof
        fdof_to_ismaster[dof] = true
      else
        ddof = -dof
        @assert (ddof_to_ddofs.ptrs[ddof+1]-ddof_to_ddofs.ptrs[ddof]) == 1
        @assert ddof_to_ddofs.data[ddof_to_ddofs.ptrs[ddof]] == ddof
        ddof_to_ismaster[ddof] = true
      end
    end
  end

  for _ddof in 1:n_ddofs
    pini = ddof_to_ddofs.ptrs[_ddof]
    pend = ddof_to_ddofs.ptrs[_ddof+1]-1
    for p in pini:pend
      dof = ddof_to_ddofs.data[p]
      if dof > 0
        @unreachable "A Dirichlet dof can only depend on Dirichlet dofs"
      else
        ddof = -dof
        @assert (ddof_to_ddofs.ptrs[ddof+1]-ddof_to_ddofs.ptrs[ddof]) == 1
        @assert ddof_to_ddofs.data[ddof_to_ddofs.ptrs[ddof]] == ddof
        ddof_to_ismaster[ddof] = true
      end
    end
  end

  fmdof_to_fdof = findall(fdof_to_ismaster)
  dmdof_to_ddof = findall(ddof_to_ismaster)

  fmdof_to_fdof, dmdof_to_ddof
end



function num_free_dofs(f::FESpaceWithLinearConstraints)
  f.num_mdofs
end

function zero_free_values(::Type{T},f::FESpaceWithLinearConstraints) where T
  zeros(T,num_free_dofs(f))
end

function get_cell_dofs(f::FESpaceWithLinearConstraints)
  f.cell_to_lmdof_to_mddof
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

function get_constraint_kernel_matrix_cols(f::FESpaceWithLinearConstraints)
  @notimplemented
end

function get_constraint_kernel_matrix_rows(f::FESpaceWithLinearConstraints)
  @notimplemented
end

function get_constraint_kernel_vector(f::FESpaceWithLinearConstraints)
  k_space = get_constraint_kernel_vector(f.space)
  k = LinearConstraintsVectorKernel(f.cell_to_lmdof_to_ldofs_and_coeffs,k_space)
end

struct LinearConstraintsVectorKernel{A,B} <: Kernel
  cell_to_lmdof_to_ldofs_and_coeffs::A
  k_space::B
end

function kernel_cache(k::LinearConstraintsVectorKernel,vec::AbstractVector,cellid::Integer)
  a = array_cache(k.cell_to_lmdof_to_ldofs_and_coeffs)
  b = kernel_cache(k.k_space,vec,cellid)
  ldof_to_v = apply_kernel!(b,k.k_space,vec,cellid)
  c = CachedArray(copy(ldof_to_v))
  a,b,c
end

function kernel_return_type(k::LinearConstraintsVectorKernel,vec::AbstractVector,cellid::Integer)
  kernel_return_type(k.k_space,vec,cellid)
end

@inline function apply_kernel!(cache,k::LinearConstraintsVectorKernel,vec::AbstractVector,cellid::Integer)
  a, b, c = cache
  ldof_to_v = apply_kernel!(b,k.k_space,vec,cellid)
  lmdof_to_ldofs_and_coeffs = getindex!(a,k.cell_to_lmdof_to_ldofs_and_coeffs,cellid)
  num_lmdofs = length(lmdof_to_ldofs_and_coeffs)
  setsize!(c,(num_lmdofs,))
  lmdof_to_v = c.array
  fill!(lmdof_to_v,zero(eltype(lmdof_to_v)))
  for lmdof in 1:num_lmdofs
    for ldof, coeff in lmdof_to_ldofs_and_coeffs[lmdof]
      lmdof_to_v[lmdof] += ldof_to_v[ldof]*coeff
    end
  end
  lmdof_to_v
end

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


