
# Background:
#
# We build a novel fe space form a given fe space plus a set of linear constraints.
# We accept any single field fe space as input.  In particular, the given fe space can also be defied
# via constraints.
#
# Assumptions:
#
#  - A constrained dof depends only on master dofs, i.e., only one level of constraints.
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
#  - sDOF a slave dof (either free or Dirichlet). It is possible to tell if it is free or Dirichlet
#    thanks to the vector sDOF_to_dof.
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

struct FESpaceWithLinearConstraints{S<:SingleFieldFESpace} <: SingleFieldFESpace
  space::S
  mDOF_to_dof::Vector
  sDOF_to_dof::Vector
  sDOF_to_mdofs::Table
  sDOF_to_coeffs::Table
  cell_to_mdofs::Table
  n_free::Int
end

# Constructors

function FESpaceWithLinearConstraints(
  space::SingleFieldFESpace,
  mDOF_to_dof::AbstractVector{<:Integer},
  sDOF_to_dof::AbstractVector{<:Integer},
  sDOF_to_mdofs::Table,
  sDOF_to_coeffs::Table,
  n_fmdofs::Int = _count_free_mdofs(mDOF_to_dof,sDOF_to_mdofs)
)
  cell_to_mdofs = _generate_cell_to_mdofs(
    space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, n_fmdofs
  )
  @check _check_constraints(space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, n_fmdofs)
  return FESpaceWithLinearConstraints(
    space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs,
    sDOF_to_coeffs, cell_to_mdofs, n_fmdofs
  )
end

function FESpaceWithLinearConstraints(
  sDOF_to_dof::AbstractVector{<:Integer},
  sDOF_to_dofs::Table,
  sDOF_to_coeffs::Table,
  space::SingleFieldFESpace
)
  mDOF_to_dof, sDOF_to_mdofs, nfmdofs = _find_master_dofs(
    sDOF_to_dof, sDOF_to_dofs, space
  )
  return FESpaceWithLinearConstraints(
    space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, nfmdofs 
  )
end

function FESpaceWithLinearConstraints(
  DOF_to_dofs::Table, DOF_to_coeffs::Table, space::SingleFieldFESpace
)
  sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = _filter_constraints(
    DOF_to_dofs, DOF_to_coeffs, space
  )
  return FESpaceWithLinearConstraints(
    sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, space
  )
end

function FESpaceWithLinearConstraints(
  fdof_to_dofs::Table,
  fdof_to_coeffs::Table,
  ddof_to_dofs::Table,
  ddof_to_coeffs::Table,
  space::SingleFieldFESpace)

  DOF_to_dofs, DOF_to_coeffs = _merge_free_and_diri_constraints(
    fdof_to_dofs, fdof_to_coeffs, ddof_to_dofs, ddof_to_coeffs
  )
  return FESpaceWithLinearConstraints(DOF_to_dofs,DOF_to_coeffs,space)
end

function _dof_to_DOF(dof,n_fdofs)
  if dof > 0
    DOF = dof
  else
    DOF = n_fdofs - dof
  end
end

function _DOF_to_dof(DOF,n_fdofs)
  if DOF > n_fdofs
    dof = -(DOF-n_fdofs)
  else
    dof = DOF
  end
end

# Implementation of FESpace interface

get_triangulation(f::FESpaceWithLinearConstraints) = get_triangulation(f.space)

get_free_dof_ids(f::FESpaceWithLinearConstraints) = Base.OneTo(f.n_free)
get_vector_type(f::FESpaceWithLinearConstraints) = get_vector_type(f.space)

get_fe_basis(f::FESpaceWithLinearConstraints) = get_fe_basis(f.space)
get_trial_fe_basis(f::FESpaceWithLinearConstraints) = get_trial_fe_basis(f.space)

CellField(f::FESpaceWithLinearConstraints,cellvals) = CellField(f.space,cellvals)

# Implementation of the SingleFieldFESpace interface

get_cell_dof_ids(f::FESpaceWithLinearConstraints) = f.cell_to_mdofs
get_fe_dof_basis(f::FESpaceWithLinearConstraints) = get_fe_dof_basis(f.space)
get_dof_value_type(f::FESpaceWithLinearConstraints) = get_dof_value_type(f.space)

get_dirichlet_dof_ids(f::FESpaceWithLinearConstraints) = Base.OneTo(length(f.mDOF_to_dof) - f.n_free)
num_dirichlet_tags(f::FESpaceWithLinearConstraints) = num_dirichlet_tags(f.space)

function get_dirichlet_dof_tag(f::FESpaceWithLinearConstraints)
  ddof_to_tag = get_dirichlet_dof_tag(f.space)
  dmdof_to_tag = zeros(eltype(ddof_to_tag),num_dirichlet_dofs(f))
  _ddof_to_dmdof_vals!(
    dmdof_to_tag,ddof_to_tag,f.mDOF_to_dof
  )
  return dmdof_to_tag
end

function get_dirichlet_dof_values(f::FESpaceWithLinearConstraints)
  ddof_to_tag = get_dirichlet_dof_values(f.space)
  dmdof_to_tag = zeros(eltype(ddof_to_tag),num_dirichlet_dofs(f))
  _ddof_to_dmdof_vals!(
    dmdof_to_tag,ddof_to_tag,f.mDOF_to_dof
  ) 
  return dmdof_to_tag
end

function scatter_free_and_dirichlet_values(f::FESpaceWithLinearConstraints,fmdof_to_val,dmdof_to_val)
  fdof_to_val = zero_free_values(f.space)
  ddof_to_val = zero_dirichlet_values(f.space)
  _mdof_to_dof_vals!(
    fdof_to_val, ddof_to_val, fmdof_to_val, dmdof_to_val,
    f.mDOF_to_dof, f.sDOF_to_dof, f.sDOF_to_mdofs, f.sDOF_to_coeffs
  )
  scatter_free_and_dirichlet_values(f.space,fdof_to_val,ddof_to_val)
end

function gather_free_and_dirichlet_values(f::FESpaceWithLinearConstraints,cell_to_ldof_to_val)
  fdof_to_val, ddof_to_val = gather_free_and_dirichlet_values(f.space,cell_to_ldof_to_val)
  fmdof_to_val = zero_free_values(f)
  dmdof_to_val = zero_dirichlet_values(f)
  _dof_to_mdof_vals!(
    fmdof_to_val, dmdof_to_val, fdof_to_val, ddof_to_val, f.mDOF_to_dof
  )
  return fmdof_to_val, dmdof_to_val
end

function gather_free_and_dirichlet_values!(fmdof_to_val,dmdof_to_val,f::FESpaceWithLinearConstraints,cell_to_ldof_to_val)
  fdof_to_val, ddof_to_val = gather_free_and_dirichlet_values(f.space,cell_to_ldof_to_val)
  _dof_to_mdof_vals!(
    fmdof_to_val, dmdof_to_val, fdof_to_val, ddof_to_val, f.mDOF_to_dof
  )
  return fmdof_to_val, dmdof_to_val
end

function _mdof_to_dof_vals!(
  fdof_to_val,ddof_to_val,fmdof_to_val,dmdof_to_val,
  mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs
)
  T = eltype(fdof_to_val)
  n_mdofs = length(mDOF_to_dof)
  n_fmdofs = length(fmdof_to_val)

  # Map free master dofs
  for (mdof,mDOF) in enumerate(1:n_fmdofs)
    dof = mDOF_to_dof[mDOF]
    if dof > 0
      fdof_to_val[dof] = fmdof_to_val[mdof]
    elseif dof < 0
      ddof_to_val[-dof] = fmdof_to_val[mdof]
    end
  end

  # Map dirichlet master dofs
  for (mdof,mDOF) in enumerate((n_fmdofs+1):n_mdofs)
    dof = mDOF_to_dof[mDOF]
    if dof < 0
      ddof_to_val[-dof] = dmdof_to_val[mdof]
    end
  end

  # Map slave dofs
  ptrs = sDOF_to_mdofs.ptrs
  for (sDOF,dof) in enumerate(sDOF_to_dof)
    val = zero(T)
    for k in ptrs[sDOF]:(ptrs[sDOF+1]-1)
      mdof = sDOF_to_mdofs.data[k]
      coeff = sDOF_to_coeffs.data[k]
      if mdof > 0
        val += fmdof_to_val[mdof]*coeff
      else
        val += dmdof_to_val[-mdof]*coeff
      end
    end
    if dof > 0
      fdof_to_val[dof] = val
    else
      ddof_to_val[-dof] = val
    end
  end

  return fdof_to_val, ddof_to_val
end

function _dof_to_mdof_vals!(
  fmdof_to_val, dmdof_to_val, fdof_to_val, ddof_to_val, mDOF_to_dof
)
  n_mdofs = length(mDOF_to_dof)
  n_fmdofs = n_mdofs - length(dmdof_to_val)

  # Map free master dofs
  for (mdof,mDOF) in enumerate(1:n_fmdofs)
    dof = mDOF_to_dof[mDOF]
    if dof > 0
      fmdof_to_val[mdof] = fdof_to_val[dof]
    elseif dof < 0
      fmdof_to_val[mdof] = ddof_to_val[-dof]
    end
  end

  # Map dirichlet master dofs
  for (mdof,mDOF) in enumerate((n_fmdofs+1):n_mdofs)
    dof = mDOF_to_dof[mDOF]
    if dof < 0
      dmdof_to_val[mdof] = ddof_to_val[-dof]
    end
  end

  return fmdof_to_val, dmdof_to_val  
end

function _ddof_to_dmdof_vals!(
  dmdof_to_val, ddof_to_val, mDOF_to_dof
)
  n_mdofs = length(mDOF_to_dof)
  n_fmdofs = n_mdofs - length(dmdof_to_val)
  for (mdof, mDOF) in enumerate((n_fmdofs+1):n_mdofs)
    dof = mDOF_to_dof[mDOF]
    if dof < 0
      dmdof_to_val[mdof] = ddof_to_val[-dof]
    end
  end
  return dmdof_to_val
end

# Constraints

ConstraintStyle(::Type{<:FESpaceWithLinearConstraints}) = Constrained()

function get_cell_isconstrained(f::FESpaceWithLinearConstraints)
  cell_dofs = get_cell_dof_ids(f.space)
  cell_mdofs = get_cell_dof_ids(f)
  mDOF_to_dof = f.mDOF_to_dof

  n_cells = length(cell_dofs)
  n_mfdofs = num_free_dofs(f)
  cell_isconstrained = Vector{Bool}(undef,n_cells)
  for cell in 1:n_cells
    dofs = view(cell_dofs,cell)
    mdofs = view(cell_mdofs,cell)

    i = 1
    mask = (length(dofs) != length(mdofs))
    while !mask && (i < length(mdofs))
      mdof = mdofs[i]
      mDOF = _dof_to_DOF(mdof,n_mfdofs)
      dof  = mDOF_to_dof[mDOF]
      mask = (dof != dofs[i])
      i += 1
    end

    cell_isconstrained[cell] = mask
  end

  return cell_isconstrained
end

function get_cell_constraints(f::FESpaceWithLinearConstraints)
  DOF_to_msDOF = generate_DOF_to_msDOF_map(f.space,f.mDOF_to_dof,f.sDOF_to_dof)
  k = LinearConstraintsMap(
    DOF_to_msDOF, f.sDOF_to_mdofs, f.sDOF_to_coeffs,
    length(f.mDOF_to_dof), num_free_dofs(f), num_free_dofs(f.space)
  )
  cell_to_mat = get_cell_constraints(f.space)
  return lazy_map(k,get_cell_dof_ids(f),get_cell_dof_ids(f.space),cell_to_mat)
end

struct LinearConstraintsMap{A,B} <: Map
  DOF_to_msDOF::Vector{Int}
  sDOF_to_mdofs::A
  sDOF_to_coeffs::B
  n_mdofs::Int
  n_fmdofs::Int
  n_fdofs::Int
end

function return_cache(k::LinearConstraintsMap,mdofs,dofs,mat)
  T = typeof(zero(eltype(mat))*zero(eltype(k.sDOF_to_coeffs.data)))
  n_lmdofs = length(mdofs)
  n_ldofs = length(dofs)
  @assert n_ldofs == size(mat,1)
  m1 = CachedArray(zeros(T,(n_lmdofs,n_ldofs)))
  m2 = CachedArray(zeros(T,(n_lmdofs,size(mat,2))))
  mdof_to_lmdof = Dict{Int,Int}()
  return m1, m2, mdof_to_lmdof
end

function evaluate!(cache,k::LinearConstraintsMap,mdofs,dofs,mat)
  m1, m2, mdof_to_lmdof = cache
  
  n_lmdofs = length(mdofs)
  n_ldofs = length(dofs)
  setsize!(m1,(n_lmdofs,n_ldofs))
  setsize!(m2,(n_lmdofs,size(mat,2)))
  a1 = get_array(m1)
  a2 = get_array(m2)
  fill!(a1,zero(eltype(a1)))

  # Precompute mdof to lmdof map
  empty!(mdof_to_lmdof)
  for (lmdof,mdof) in enumerate(mdofs)
    mdof_to_lmdof[mdof] = lmdof
  end

  # Precompute DOF to msDOF map
  o = one(eltype(a1))
  ptrs = k.sDOF_to_mdofs.ptrs
  for (ldof,dof) in enumerate(dofs)
    DOF = _dof_to_DOF(dof,k.n_fdofs)
    msDOF = k.DOF_to_msDOF[DOF]
    if msDOF <= k.n_mdofs # master dof
      mDOF = msDOF
      mdof = _DOF_to_dof(mDOF,k.n_fmdofs)
      lmdof = mdof_to_lmdof[mdof]
      a1[lmdof,ldof] = o
    else # slave dof
      sDOF = msDOF - k.n_mdofs
      for i in ptrs[sDOF]:(ptrs[sDOF+1]-1)
        mdof  = k.sDOF_to_mdofs.data[i]
        coeff = k.sDOF_to_coeffs.data[i]
        lmdof = mdof_to_lmdof[mdof]
        a1[lmdof,ldof] = coeff
      end
    end
  end

  #TODO this is not always needed
  mul!(a2,a1,mat)
  a2
end

# Private methods

function _check_constraints(
  space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, n_fmdofs
)
  n_fdofs = num_free_dofs(space)
  n_dofs = n_fdofs + num_dirichlet_dofs(space)
  DOF_is_master = fill(false,n_dofs)
  for dof in mDOF_to_dof
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_is_master[DOF] = true
  end

  mDOF_is_master = fill(true,length(mDOF_to_dof))
  for (mDOF,dof) in enumerate(mDOF_to_dof)
    DOF = _dof_to_DOF(dof,n_fdofs)
    mDOF_is_master[mDOF] = DOF_is_master[DOF]
  end

  c = array_cache(sDOF_to_mdofs)
  for (sDOF,dof) in enumerate(sDOF_to_dof)
    DOF = _dof_to_DOF(dof,n_fdofs)
    @check !DOF_is_master[DOF] "A dof cannot be both master and slave."
    mdofs = getindex!(c,sDOF_to_mdofs,sDOF)
    for mdof in mdofs
      mDOF = _dof_to_DOF(mdof,n_fmdofs)
      @check mDOF_is_master[mDOF] "Constraints cannot be recursive."
      @check (dof > 0) || (mdof < 0) "Dirichlet dofs can only be constrained by Dirichlet dofs."
    end
  end

  return true
end

function _filter_constraints(DOF_to_dofs, DOF_to_coeffs, space)
  n_fdofs = num_free_dofs(space)
  isslave(DOF,dofs) = (dofs != [_DOF_to_dof(DOF,n_fdofs)])  
  sDOF_to_DOF = findall(map(isslave,eachindex(DOF_to_dofs),DOF_to_dofs))
  sDOF_to_dof = map(Base.Fix2(_DOF_to_dof,n_fdofs),sDOF_to_DOF)
  sDOF_to_dofs = DOF_to_dofs[sDOF_to_DOF]
  sDOF_to_coeffs = DOF_to_coeffs[sDOF_to_DOF]
  return sDOF_to_dof, sDOF_to_dofs, Table(sDOF_to_coeffs.data,sDOF_to_dofs.ptrs)
end

function _merge_free_and_diri_constraints(fdof_to_dofs, fdof_to_coeffs, ddof_to_dofs, ddof_to_coeffs)
  DOF_to_dofs = append_tables_globally(fdof_to_dofs,ddof_to_dofs)
  DOF_to_coeffs = Table(
    vcat(fdof_to_coeffs.data,ddof_to_coeffs.data), DOF_to_dofs.ptrs
  )
  return DOF_to_dofs, DOF_to_coeffs
end

function _count_free_mdofs(mDOF_to_dof,sDOF_to_mdofs)
  a = count(>(0), mDOF_to_dof; init=0)
  b = maximum(sDOF_to_mdofs.data)
  return max(a, b)
end

function _find_master_dofs(
  sDOF_to_dof, sDOF_to_dofs, space
)
  n_fdofs = num_free_dofs(space)
  n_dofs = n_fdofs + num_dirichlet_dofs(space)

  # Flag master dofs 
  DOF_ismaster = fill(true,n_dofs)
  for dof in sDOF_to_dof
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_ismaster[DOF] = false
  end

  # Counting mdofs
  n_fmdofs = count(view(DOF_ismaster,1:n_fdofs))
  n_dmdofs = count(view(DOF_ismaster,(n_fdofs+1):n_dofs))
  n_mdofs = n_fmdofs + n_dmdofs

  # DOF to mdof mapping
  kf, kd = 1, 1
  mDOF_to_dof = Vector{Int32}(undef,n_mdofs)
  DOF_to_mdof = Vector{Int32}(undef,n_dofs)
  for DOF in 1:n_dofs
    if DOF_ismaster[DOF]
      dof = _DOF_to_dof(DOF,n_fdofs)
      if dof > 0
        mDOF_to_dof[kf] = dof
        DOF_to_mdof[DOF] = kf
        kf += 1
      else
        mDOF_to_dof[n_fmdofs + kd] = dof
        DOF_to_mdof[DOF] = -kd
        kd += 1
      end
    end
  end

  # Map constraints
  data = zeros(Int32,length(sDOF_to_dofs.data))
  for (i,dof) in enumerate(sDOF_to_dofs.data)
    DOF = _dof_to_DOF(dof,n_fdofs)
    data[i] = DOF_to_mdof[DOF]
  end
  sDOF_to_mdofs = Table(data,sDOF_to_dofs.ptrs)

  return mDOF_to_dof, sDOF_to_mdofs, n_fmdofs
end

function _generate_cell_to_mdofs(
  space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, n_fmdofs
)
  DOF_to_msDOF = generate_DOF_to_msDOF_map(
    space,mDOF_to_dof,sDOF_to_dof
  )
  cell_dofs = get_cell_dof_ids(space)
  c1 = array_cache(cell_dofs)
  c2 = array_cache(sDOF_to_mdofs)

  n_cells = length(cell_dofs)
  n_fdofs = num_free_dofs(space)
  n_mDOFs = length(mDOF_to_dof)

  acc = OrderedSet{Int32}()
  ptrs = zeros(Int32,n_cells+1)
  for cell in 1:n_cells
    dofs = getindex!(c1,cell_dofs,cell)
    for dof in dofs
      DOF = _dof_to_DOF(dof,n_fdofs)
      msDOF = DOF_to_msDOF[DOF]
      if msDOF <= n_mDOFs # master dof
        mdof = _DOF_to_dof(msDOF,n_fmdofs)
        push!(acc, mdof)
      else # slave dof
        sDOF = msDOF - n_mDOFs
        mdofs = getindex!(c2,sDOF_to_mdofs,sDOF)
        push!(acc, mdofs...)
      end
    end
    ptrs[cell+1] += length(acc)
    empty!(acc)
  end
  Arrays.length_to_ptrs!(ptrs)

  data = zeros(Int32,ptrs[end]-1)
  for cell in 1:n_cells
    dofs = getindex!(c1,cell_dofs,cell)
    for dof in dofs
      DOF = _dof_to_DOF(dof,n_fdofs)
      msDOF = DOF_to_msDOF[DOF]
      if msDOF <= n_mDOFs # master dof
        mdof = _DOF_to_dof(msDOF,n_fmdofs)
        push!(acc, mdof)
      else # slave dof
        sDOF = msDOF - n_mDOFs
        mdofs = getindex!(c2,sDOF_to_mdofs,sDOF)
        push!(acc, mdofs...)
      end
    end
    data[ptrs[cell]:(ptrs[cell+1]-1)] .= collect(acc)
    empty!(acc)
  end

  return Table(data,ptrs)
end

function generate_DOF_to_msDOF_map(space, mDOF_to_dof, sDOF_to_dof)
  n_mDOFs = length(mDOF_to_dof)
  n_fdofs = num_free_dofs(space)
  n_dofs  = n_fdofs + num_dirichlet_dofs(space)
  DOF_to_msDOF = Vector{Int}(undef,n_dofs)

  for (mDOF,dof) in enumerate(mDOF_to_dof)
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_to_msDOF[DOF] = mDOF
  end

  for (sDOF,dof) in enumerate(sDOF_to_dof)
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_to_msDOF[DOF] = sDOF + n_mDOFs
  end

  return DOF_to_msDOF
end
