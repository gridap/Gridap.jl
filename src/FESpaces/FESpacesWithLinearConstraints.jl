
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
  n_fdofs = length(fdof_to_dofs)
  DOF_to_DOFs = append_tables_globally(fdof_to_dofs,ddof_to_dofs)
  for i in 1:length(DOF_to_DOFs.data)
    dof = DOF_to_DOFs.data[i]
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_to_DOFs.data[i] = DOF
  end
  DOF_to_coeffs = Table(vcat(fdof_to_coeffs.data,ddof_to_coeffs.data),DOF_to_DOFs.ptrs)
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
  cell_to_lmdof_to_mdof_ptrs = zeros(eltype(cell_to_ldof_to_dof.ptrs),n_cells+1)

  for cell in 1:n_cells
    mdofs = Set{Int}()
    pini = cell_to_ldof_to_dof.ptrs[cell]
    pend = cell_to_ldof_to_dof.ptrs[cell+1]-1
    for p in pini:pend
      dof = cell_to_ldof_to_dof.data[p]
      DOF = _dof_to_DOF(dof,n_fdofs)
      qini = DOF_to_mDOFs.ptrs[DOF]
      qend = DOF_to_mDOFs.ptrs[DOF+1]-1
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
      qini = DOF_to_mDOFs.ptrs[DOF]
      qend = DOF_to_mDOFs.ptrs[DOF+1]-1
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

function get_dirichlet_dof_tag(f::FESpaceWithLinearConstraints)
  ddof_to_tag = get_dirichlet_dof_tag(f.space)
  dmdof_to_tag = zeros(eltype(ddof_to_tag),num_dirichlet_dofs(f))
  _setup_ddof_to_tag!(
    dmdof_to_tag,
    ddof_to_tag,
    f.mDOF_to_DOF,
    f.n_fdofs,
    f.n_fmdofs)
  dmdof_to_tag
end

function _setup_ddof_to_tag!(
  dmdof_to_tag,
  ddof_to_tag,
  mDOF_to_DOF,
  n_fdofs,
  n_fmdofs)

  for mDOF in (n_fmdofs+1):length(mDOF_to_DOF)
    mdof = _DOF_to_dof(mDOF,n_fmdofs)
    @assert mdof < 0 "Dirichlet dofs can only depend on Dirichlet dofs"
    dmdof = -mdof
    DOF = mDOF_to_DOF[mDOF]
    dof = _DOF_to_dof(DOF,n_fdofs)
    @assert dof < 0 "Dirichlet dofs can only depend on Dirichlet dofs"
    ddof = -dof
    dmdof_to_tag[dmdof] = ddof_to_tag[ddof]
  end
end

function get_dirichlet_values(f::FESpaceWithLinearConstraints)
  ddof_to_tag = get_dirichlet_values(f.space)
  dmdof_to_tag = zeros(eltype(ddof_to_tag),num_dirichlet_dofs(f))
  _setup_ddof_to_tag!(
    dmdof_to_tag,
    ddof_to_tag,
    f.mDOF_to_DOF,
    f.n_fdofs,
    f.n_fmdofs)
  dmdof_to_tag
end

function scatter_free_and_dirichlet_values(f::FESpaceWithLinearConstraints,fmdof_to_val,dmdof_to_val)
  fdof_to_val = zero_free_values(eltype(fmdof_to_val),f.space)
  ddof_to_val = zero_dirichlet_values(f.space)
  _setup_dof_to_val!(
    fdof_to_val,
    ddof_to_val,
    fmdof_to_val,
    dmdof_to_val,
    f.DOF_to_mDOFs,
    f.DOF_to_coeffs,
    f.n_fdofs,
    f.n_fmdofs)
  scatter_free_and_dirichlet_values(f.space,fdof_to_val,ddof_to_val)
end

function _setup_dof_to_val!(
  fdof_to_val,
  ddof_to_val,
  fmdof_to_val,
  dmdof_to_val,
  DOF_to_mDOFs,
  DOF_to_coeffs,
  n_fdofs,
  n_fmdofs)

  T = eltype(fdof_to_val)

  for DOF in 1:length(DOF_to_mDOFs)
    pini = DOF_to_mDOFs.ptrs[DOF]
    pend = DOF_to_mDOFs.ptrs[DOF+1]-1
    val = zero(T)
    for p in pini:pend
      mDOF = DOF_to_mDOFs.data[p]
      coeff = DOF_to_coeffs.data[p]
      mdof = _DOF_to_dof(mDOF,n_fmdofs)
      if mdof > 0
        fmdof = mdof
        val += fmdof_to_val[fmdof]*coeff
      else
        dmdof = -mdof
        val += dmdof_to_val[dmdof]*coeff
      end
    end
    dof = _DOF_to_dof(DOF,n_fdofs)
    if dof > 0
      fdof = dof
      fdof_to_val[fdof] = val
    else
      ddof = -dof
      ddof_to_val[ddof] = val
    end
  end

end

function gather_free_and_dirichlet_values(f::FESpaceWithLinearConstraints,cell_to_ludof_to_val)
  fdof_to_val, ddof_to_val = gather_free_and_dirichlet_values(f.space,cell_to_ludof_to_val)
  fmdof_to_val = zero_free_values(eltype(fdof_to_val),f)
  dmdof_to_val = zero_dirichlet_values(f)
  _setup_mdof_to_val!(
    fmdof_to_val,
    dmdof_to_val,
    fdof_to_val,
    ddof_to_val,
    f.mDOF_to_DOF,
    f.n_fdofs,
    f.n_fmdofs)
  fmdof_to_val, dmdof_to_val
end

function _setup_mdof_to_val!(
  fmdof_to_val,
  dmdof_to_val,
  fdof_to_val,
  ddof_to_val,
  mDOF_to_DOF,
  n_fdofs,
  n_fmdofs)

  for mDOF in 1:length(mDOF_to_DOF)
    DOF = mDOF_to_DOF[mDOF]
    dof = _DOF_to_dof(DOF,n_fdofs)
    if dof > 0
      fdof = dof
      val = fdof_to_val[fdof]
    else
      ddof = -dof
      val = ddof_to_val[ddof]
    end
    mdof = _DOF_to_dof(mDOF,n_fmdofs)
    if mdof > 0
      fmdof = mdof
      fmdof_to_val[fmdof] = val
    else
      dmdof = -mdof
      dmdof_to_val[dmdof] = val
    end
  end

end

# Implementation of FESpace interface

num_free_dofs(f::FESpaceWithLinearConstraints) = f.n_fmdofs

function get_cell_basis(f::FESpaceWithLinearConstraints)
  get_cell_basis(f.space)
end

function zero_free_values(::Type{T},f::FESpaceWithLinearConstraints) where T
  zeros(T,num_free_dofs(f))
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
  LinearConstraintsVectorKernel(
    f.cell_to_lmdof_to_mdof,
    f.cell_to_ldof_to_dof,
    f.DOF_to_mDOFs,
    f.DOF_to_coeffs,
    k_space,
    length(f.mDOF_to_DOF),
    f.n_fmdofs,
    f.n_fdofs)
end

struct LinearConstraintsVectorKernel{A,B,C} <: Kernel
  cell_to_lmdof_to_mdof::A
  cell_to_ldof_to_dof::A
  DOF_to_mDOFs::A
  DOF_to_coeffs::B
  k_space::C
  n_mDOFs::Int
  n_fmdofs::Int
  n_fdofs::Int
end

function kernel_cache(k::LinearConstraintsVectorKernel,ludof_to_val::AbstractVector,cellid::Integer)
  a = kernel_cache(k.k_space, ludof_to_val, cellid)
  ldof_to_val = apply_kernel!(a, k.k_space, ludof_to_val, cellid)
  b = CachedArray(copy(ldof_to_val))
  mDOF_to_lmdof = zeros(Int8,length(k.n_mDOFs))
  a, b, mDOF_to_lmdof
end

function kernel_return_type(k::LinearConstraintsVectorKernel,ludof_to_val::AbstractVector,cellid::Integer)
  kernel_return_type(k.k_space,ludof_to_val,cellid)
end

@inline function apply_kernel!(
  cache,k::LinearConstraintsVectorKernel,ludof_to_val::AbstractVector,cellid::Integer)

  a, b, mDOF_to_lmdof = cache

  ldof_to_val = apply_kernel!(a, k.k_space, ludof_to_val, cellid)

  pini = k.cell_to_lmdof_to_mdof.ptrs[cellid]
  pend = k.cell_to_lmdof_to_mdof.ptrs[cellid+1]-1
  for (lmdof,p) in enumerate(pini:pend)
    mdof = k.cell_to_lmdof_to_mdof.data[p]
    mDOF = _dof_to_DOF(mdof,k.n_fmdofs)
    mDOF_to_lmdof[mDOF] = lmdof
  end

  n_lmdofs = length(pini:pend)
  setsize!(b,(n_lmdofs,))
  lmdof_to_val = b.array
  fill!(lmdof_to_val,zero(eltype(lmdof_to_val)))

  pini = k.cell_to_ldof_to_dof.ptrs[cellid]
  pend = k.cell_to_ldof_to_dof.ptrs[cellid+1]-1
  for (ldof,p) in enumerate(pini:pend)
    val = ldof_to_val[ldof]
    dof = k.cell_to_ldof_to_dof.data[p]
    DOF = _dof_to_DOF(dof,k.n_fdofs)
    qini = k.DOF_to_mDOFs.ptrs[DOF]
    qend = k.DOF_to_mDOFs.ptrs[DOF+1]-1
    for q in qini:qend
      mDOF = k.DOF_to_mDOFs.data[q]
      coeff = k.DOF_to_coeffs.data[q]
      lmdof = mDOF_to_lmdof[mDOF]
      lmdof_to_val[lmdof] += coeff*val
    end
  end

  lmdof_to_val
end

