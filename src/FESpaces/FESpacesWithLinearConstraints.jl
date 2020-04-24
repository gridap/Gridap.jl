
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
#   The sign of the value tells if it is one or the other.
#
#  - DOF: a dof (either free or Dirichlet) in the original space.
#    In this case, Dirichlet dofs are also represented with a positive id "pas the end" of free dof ids.
#    Useful for indexing arrays.
#    It will refer to a free dof if 0 < DOF <= n_fdofs, and Dirichlet if n_fdofs < DOF <= n_dofs
#    with n_ddofs=n_dofs-n_fdofs
#
#  - fdof: a free dof in the original space. We have: fdof = dof = DOF
#
#  - ddof: a Dirichlet dof in the original space. We have: ddof = -dof = DOF - n_fdofs
#
#  - mdof: a master dof (i.e. a dof in the constrained space).
#    It can be free or Dirichlet. The sign of the value tells if it is one or the other.
#
#  - mDOF: a master dof (i.e. a dof in the constrained space).
#    In this case, Dirichlet dofs are also represented with a positive id "pas the end" of free dof ids.
#    It will refer to a free dof if 0 < mDOF <= n_fmdofs, and Dirichlet if n_fmdofs < mDOF <= n_mdofs
#    with n_dmdofs=n_mdofs-n_fmdofs
#
#  - fmdof: a free master dof. We have fmdof = mdof = mDOF
#
#  - dmdof: a Dirichlet master dof. We have dmdof = -mdof = mDOF - n_fmdofs
#
#  - ldof: a local dof in a cell for the original space. Always positive.
#
#  - lmdof: local master dof in a cell. That is, a local dofs in a cell for the constrained space.
#    Always positive.
#
#  - ludof: A local dof in the unconstrained interpolation space of the underlying cell. Always positive.
#
#  - coeff: A coefficient in a linear constraint

struct FESpaceWithLinearConstraints <: SingleFieldFESpace
  space::SingleFieldFESpace
  n_fdofs::Int
  n_fmdofs::Int
  mDOF_to_DOF::Vector
  DOF_to_mDOFs::Table
  DOF_to_coeffs::Table
  cell_to_lmdof_to_mdof::Table
  cell_to_ldof_to_dof::Table
end

# Constructors

function FESpaceWithLinearConstraints(
  fdof_to_dofs::Table,
  fdof_to_coeffs::Table,
  ddof_to_dofs::Table,
  ddof_to_coeffs::Table,
  space::SingleFieldFESpace)

  DOF_to_DOFs, DOF_to_coeffs = _merge_free_and_diri_constraints(
    fdof_to_dofs, fdof_to_coeffs, ddof_to_dofs, ddof_to_coeffs)
  FESpaceWithLinearConstraints!(DOF_to_DOFs,DOF_to_coeffs,space)
end

function _merge_free_and_diri_constraints(fdof_to_dofs, fdof_to_coeffs, ddof_to_dofs, ddof_to_coeffs)
  DOF_to_DOFs = append_tables_globally(fdof_to_dofs,ddof_to_dofs)
  DOF_to_coeffs = Table(vcat(fdof_to_coeffs,ddof_to_coeffs),DOF_to_DOFs.ptrs)
  DOF_to_DOFs, DOF_to_coeffs
end

function FESpaceWithLinearConstraints(DOF_to_DOFs::Table,DOF_to_coeffs::Table,space::SingleFieldFESpace)
  FESpaceWithLinearConstraints!(copy(DOF_to_DOFs),copy(DOF_to_coeffs),space::SingleFieldFESpace)
end

function FESpaceWithLinearConstraints!(DOF_to_DOFs::Table,DOF_to_coeffs::Table,space::SingleFieldFESpace)

  n_fdofs = num_free_dofs(space)
  mDOF_to_DOF, n_fmdofs = _find_master_dofs(DOF_to_DOFs,n_fdofs)
  DOF_to_mDOFs = _renumber_constraints!(DOF_to_DOFs,mDOF_to_DOF)
  cell_to_ldof_to_dof = Table(get_cell_dofs(space))
  cell_to_lmdof_to_mdof = _setup_cell_to_lmdof_to_mdof(cell_to_ldof_to_dof,DOF_to_mDOFs,n_fdofs,n_fmdofs)

  FESpaceWithLinearConstraints(
    space,
    n_fdofs,
    n_fmdofs,
    mDOF_to_DOF,
    DOF_to_mDOFs,
    DOF_to_coeffs,
    cell_to_lmdof_to_mdof,
    cell_to_ldof_to_dof)

end

function _find_master_dofs(DOF_to_DOFs,n_fdofs)
  n_DOFs = length(DOF_to_DOFs)
  DOF_to_ismaster = fill(false,n_DOFs)
  for DOF in 1:n_DOFs
    pini = DOF_to_DOFs.ptrs[DOF]
    pend = DOF_to_DOFs.ptrs[DOF+1]-1
    for p in pini:pend
      _DOF = DOF_to_DOFs.data[p]
      @assert (DOF_to_DOFs.ptrs[_DOF+1]-DOF_to_DOFs.ptrs[_DOF]) == 1 "Rcursive constraints not allowed now"
      @assert DOF_to_DOFs.data[DOF_to_DOFs.ptrs[_DOF]] == _DOF "Rcursive constraints not allowed now"
      DOF_to_ismaster[_DOF] = true
  end
  n_fmdofs = 0
  for DOF in 1:n_fdofs
    if DOF_to_ismaster[DOF]
      n_fmdofs += 1
    end
  end
  mDOF_to_DOF = findall(DOF_to_ismaster)
  mDOF_to_DOF, n_fmdofs
end

function _renumber_constraints!(DOF_to_DOFs,mDOF_to_DOF)
  DOF_to_mDOF = zeros(eltype(DOF_to_DOFs.data),length(DOF_to_DOFs))
  DOF_to_mDOF[mDOF_to_DOF] .= 1:length(mDOF_to_DOF)
  for i in 1:length(DOF_to_DOFs.data)
    DOF = DOF_to_DOFs.data[i]
    mDOF = DOF_to_mDOF[DOF]
    DOF_to_DOFs.data[i] = mDOF
  end
  DOF_to_DOFs
end

function _setup_cell_to_lmdof_to_mdof(cell_to_ldof_to_dof,DOF_to_mDOFs,n_fdofs,n_fmdofs)

  n_cells = length(cell_to_ldof_to_dof)
  cell_to_lmdof_to_mdof_ptrs = zeros(eltype(cell_to_ldof_to_dof.ptrs),n_cells)

  for cell in 1:n_cells
    mdofs = Set{Int}()
    pini = cell_to_ldof_to_dof.ptrs[cell]
    pend = cell_to_ldof_to_dof.ptrs[cell+1]-1
    for p in pini:pend
      dof = cell_to_ldof_to_dof.data[p]
      DOF = _dof_to_DOF(dof,n_fdofs)
      qini = DOF_to_mDOFs.ptrs[fdof]
      qend = DOF_to_mDOFs.ptrs[fdof+1]-1
      for q in qini:qend
        mDOF = DOF_to_mDOFs.data[q]
        mdof = _DOF_to_dof(mDOF,n_fmdofs)
        push!(mdofs,mdof)
      end
    end
    cell_to_lmdof_to_mdof_ptrs[cell+1] = length(mdofs)
  end

  length_to_ptrs!(cell_to_lmdof_to_mdof_ptrs)
  ndata = cell_to_lmdof_to_mdof_ptrs[end]-1
  cell_to_lmdof_to_mdof_data = zeros(eltype(cell_to_ldof_to_dof.data),ndata)

  for cell in 1:n_cells
    mdofs = Set{Int}()
    pini = cell_to_ldof_to_dof.ptrs[cell]
    pend = cell_to_ldof_to_dof.ptrs[cell+1]-1
    for p in pini:pend
      dof = cell_to_ldof_to_dof.data[p]
      DOF = _dof_to_DOF(dof,n_fdofs)
      qini = DOF_to_mDOFs.ptrs[fdof]
      qend = DOF_to_mDOFs.ptrs[fdof+1]-1
      for q in qini:qend
        mDOF = DOF_to_mDOFs.data[q]
        mdof = _DOF_to_dof(mDOF,n_fmdofs)
        push!(mdofs,mdof)
      end
    end
    o = cell_to_lmdof_to_mdof_ptrs[cell]-1
    for (lmdof, mdof) in enumerate(mdofs)
      cell_to_lmdof_to_mdof_data[o+lmdof] = mdof
    end
  end

  Table(cell_to_lmdof_to_mdof_data,cell_to_lmdof_to_mdof_ptrs)
end

@inline function _dof_to_DOF(dof,n_fdofs)
  if dof > 0
    DOF = dof
  else
    DOF = n_fdofs - dof
  end
end

@inline function _DOF_to_dof(DOF,n_fdofs)
  if DOF > n_fdofs
    dof = -(DOF-n_fdofs)
  else
    dof = DOF
  end
end

# Implementation of the SingleFieldFESpace interface

function get_cell_dofs(f::FESpaceWithLinearConstraints)
  f.cell_to_lmdof_to_mddof
end

function get_cell_dof_basis(f::FESpaceWithLinearConstraints)
  get_cell_dof_basis(f.space)
end

num_dirichlet_dofs(f::FESpaceWithLinearConstraints) = length(f.mDOF_to_DOF) - f.n_fmdofs

function zero_dirichlet_values(f::FESpaceWithLinearConstraints)
  # TODO
  zeros(Float64,num_dirichlet_dofs(f))
end

num_dirichlet_tags(f::FESpaceWithLinearConstraints) = num_dirichlet_tags(f.space)

get_dirichlet_dof_tag(f::FESpaceWithLinearConstraints) = get_dirichlet_dof_tag(f.space)

function scatter_free_and_dirichlet_values(f::FESpaceWithLinearConstraints,fmdof_to_val,dmdof_to_val)
  fdof_to_val = zero_free_values(eltype(fmdof_to_val),f.space)
  ddof_to_val = zero_dirichlet_values(f.space)
  _setup_dof_to_val_and_mdof_to_val!(fdof_to_val,mdof_to_val,ddof_to_val,fdof_to_mddofs,fdof_to_coeffs)
  scatter_free_and_dirichlet_values(f.space,fdof_to_val,ddof_to_val)
end

function _setup_dof_to_val_and_mdof_to_val!(
  fdof_to_val,ddof_to_val,DOF_to_mDOFs,n_fmdofs)



end








num_free_dofs(f::FESpaceWithLinearConstraints) = f.n_fmdofs


function zero_free_values(::Type{T},f::FESpaceWithLinearConstraints) where T
  zeros(T,num_free_dofs(f))
end








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

@inline function apply_kernel!(cache,k::LinearConstraintsVectorKernel,vec::AbstractVector,cellid::Integer)
  # prepare the map fmdof_to_lmdof and dmdof_to_lmdof
  a = k.cell_ldof_to_dofs.ptrs[cellid]-1
  for ldof in 1:n_ldofs
    dof = k.cell_ldof_to_dofs.data[a+ldof]
    if dof > 0
      fdof = dof
      pini = k.fdof_to_mdofs.ptrs[fdof]
      pend = k.fdof_to_mdofs.ptrs[fdof+1]-1
      for p in pini:pend
        mdof = k.fdof_to_mdofs.data[p]
        if mddof>0
        else
        end
      end
    else
      ddof = -dof
    end

  end
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


