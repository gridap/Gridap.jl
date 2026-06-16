
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

"""
    struct FESpaceWithLinearConstraints{S<:SingleFieldFESpace} <: SingleFieldFESpace
"""
struct FESpaceWithLinearConstraints{S<:SingleFieldFESpace} <: SingleFieldFESpace
  space::S
  mDOF_to_dof::Vector
  sDOF_to_dof::Vector
  sDOF_to_mdofs::Table
  sDOF_to_coeffs::Table
  cell_to_mdofs::Table
  dmdof_to_offsets::Vector
  n_free::Int
end

# Constructors

function FESpaceWithLinearConstraints(
  space::SingleFieldFESpace,
  mDOF_to_dof::AbstractVector{<:Integer},
  sDOF_to_dof::AbstractVector{<:Integer},
  sDOF_to_mdofs::Table,
  sDOF_to_coeffs::Table,
  n_fmdofs::Int = _count_free_mdofs(mDOF_to_dof,sDOF_to_mdofs),
  dmdof_to_offsets::AbstractVector = zeros(Float64, length(mDOF_to_dof) - n_fmdofs)
)
  cell_to_mdofs = _generate_cell_to_mdofs(
    space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, n_fmdofs
  )
  @check _check_constraints(space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, n_fmdofs)
  return FESpaceWithLinearConstraints(
    space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs,
    sDOF_to_coeffs, cell_to_mdofs, dmdof_to_offsets, n_fmdofs
  )
end

"""
    FESpaceWithLinearConstraints(
        sDOF_to_dof::AbstractVector{<:Integer},
        sDOF_to_dofs::Table,
        sDOF_to_coeffs::Table,
        space::SingleFieldFESpace
    )
"""
function FESpaceWithLinearConstraints(
  sDOF_to_dof::AbstractVector{<:Integer},
  sDOF_to_dofs::Table,
  sDOF_to_coeffs::Table,
  space::SingleFieldFESpace;
  sDOF_to_offsets=nothing
)
  mDOF_to_dof, sDOF_to_mdofs, n_fmdofs, n_dmdofs =
    _find_master_dofs(sDOF_to_dof, sDOF_to_dofs, space)
  if !isnothing(sDOF_to_offsets)
    mDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, dmdof_to_offsets =
      _attach_offsets(mDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, sDOF_to_offsets, n_dmdofs)
  else
    dmdof_to_offsets = eltype(sDOF_to_coeffs.data)[]
  end
  return FESpaceWithLinearConstraints(
    space, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, n_fmdofs, dmdof_to_offsets
  )
end

function FESpaceWithLinearConstraints(
  DOF_to_dofs::Table, DOF_to_coeffs::Table, space::SingleFieldFESpace,
  DOF_to_offsets=nothing
)
  sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets = _filter_constraints(
    DOF_to_dofs, DOF_to_coeffs, DOF_to_offsets, space
  )
  return FESpaceWithLinearConstraints(
    sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, space; sDOF_to_offsets
  )
end

function FESpaceWithLinearConstraints(
  fdof_to_dofs::Table,
  fdof_to_coeffs::Table,
  ddof_to_dofs::Table,
  ddof_to_coeffs::Table,
  space::SingleFieldFESpace
)
  DOF_to_dofs, DOF_to_coeffs = _merge_free_and_diri_constraints(
    fdof_to_dofs, fdof_to_coeffs, ddof_to_dofs, ddof_to_coeffs
  )
  return FESpaceWithLinearConstraints(DOF_to_dofs, DOF_to_coeffs, space)
end

_dof_to_DOF(dof, n_fdofs) = ifelse(iszero(dof), 0, ifelse(dof > 0, dof, n_fdofs - dof))
_DOF_to_dof(DOF, n_fdofs) = ifelse(DOF > n_fdofs, -(DOF - n_fdofs), DOF)

function _dof_to_DOF!(dofs::Vector{Int}, n_fdofs::Int)
  @inbounds for i in eachindex(dofs)
    dofs[i] = _dof_to_DOF(dofs[i], n_fdofs)
  end
  return dofs
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
  ddof_to_val = get_dirichlet_dof_values(f.space)
  dmdof_to_val = zeros(eltype(ddof_to_val), num_dirichlet_dofs(f))
  _ddof_to_dmdof_vals!(dmdof_to_val, ddof_to_val, f.mDOF_to_dof)
  n_offsets = length(f.dmdof_to_offsets)
  dmdof_to_val[end-n_offsets+1:end] .= f.dmdof_to_offsets
  return dmdof_to_val
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
    while !mask && (i <= length(mdofs))
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

function LinearConstraintsMap(
  sDOF_to_dof::AbstractVector{<:Integer},
  sDOF_to_dofs::Table,
  sDOF_to_coeffs::Table,
  space::SingleFieldFESpace
)
  mDOF_to_dof, sDOF_to_mdofs, n_fmdofs, _ =
    _find_master_dofs(sDOF_to_dof, sDOF_to_dofs, space)
  DOF_to_msDOF = generate_DOF_to_msDOF_map(space,mDOF_to_dof,sDOF_to_dof)
  return LinearConstraintsMap(
    DOF_to_msDOF, sDOF_to_mdofs, sDOF_to_coeffs,
    length(mDOF_to_dof), n_fmdofs, num_free_dofs(space)
  )
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
    iszero(dof) && continue
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_is_master[DOF] = true
  end

  mDOF_is_master = fill(true,length(mDOF_to_dof))
  for (mDOF,dof) in enumerate(mDOF_to_dof)
    DOF = _dof_to_DOF(dof,n_fdofs)
    mDOF_is_master[mDOF] = iszero(DOF) || DOF_is_master[DOF]
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

function _filter_constraints(DOF_to_dofs, DOF_to_coeffs, DOF_to_offsets, space)
  n_fdofs = num_free_dofs(space)
  isslave(DOF,dofs) = (dofs != [_DOF_to_dof(DOF,n_fdofs)])
  sDOF_to_DOF = findall(map(isslave,eachindex(DOF_to_dofs),DOF_to_dofs))
  sDOF_to_dof = map(Base.Fix2(_DOF_to_dof,n_fdofs),sDOF_to_DOF)
  sDOF_to_dofs    = DOF_to_dofs[sDOF_to_DOF]
  sDOF_to_coeffs  = DOF_to_coeffs[sDOF_to_DOF]
  sDOF_to_offsets = isnothing(DOF_to_offsets) ? nothing : DOF_to_offsets[sDOF_to_DOF]
  return sDOF_to_dof, sDOF_to_dofs, Table(sDOF_to_coeffs.data,sDOF_to_dofs.ptrs), sDOF_to_offsets
end

function _merge_free_and_diri_constraints(
  fdof_to_dofs, fdof_to_coeffs, ddof_to_dofs, ddof_to_coeffs
)
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

function _find_master_dofs(sDOF_to_dof, sDOF_to_dofs, space)
  n_fdofs = num_free_dofs(space)
  n_dofs  = n_fdofs + num_dirichlet_dofs(space)

  DOF_ismaster = fill(true, n_dofs)
  for dof in sDOF_to_dof
    DOF = _dof_to_DOF(dof, n_fdofs)
    DOF_ismaster[DOF] = false
  end

  n_fmdofs = count(view(DOF_ismaster, 1:n_fdofs))
  n_dmdofs = count(view(DOF_ismaster, (n_fdofs+1):n_dofs))

  kf, kd = 1, 1
  mDOF_to_dof = Vector{Int32}(undef, n_fmdofs + n_dmdofs)
  DOF_to_mdof = Vector{Int32}(undef, n_dofs)
  for DOF in 1:n_dofs
    if DOF_ismaster[DOF]
      dof = _DOF_to_dof(DOF, n_fdofs)
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

  data = zeros(Int32, length(sDOF_to_dofs.data))
  for (i, dof) in enumerate(sDOF_to_dofs.data)
    data[i] = DOF_to_mdof[_dof_to_DOF(dof, n_fdofs)]
  end

  return mDOF_to_dof, Table(data, sDOF_to_dofs.ptrs), n_fmdofs, n_dmdofs
end

# Extends mDOF_to_dof, sDOF_to_mdofs, and sDOF_to_coeffs with fictitious Dirichlet
# masters encoding affine inhomogeneities. One fictitious master per nonzero offset,
# appended with coefficient 1. Returns (mDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, dmdof_to_offsets).
function _attach_offsets(mDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, sDOF_to_offsets, n_dmdofs)
  n_slaves  = length(sDOF_to_mdofs)
  n_offsets = count(!iszero, sDOF_to_offsets)
  OT = eltype(sDOF_to_offsets)
  CT = eltype(sDOF_to_coeffs.data)
  iszero(n_offsets) && return mDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, OT[]

  dmdof_to_offsets = Vector{OT}(undef, n_offsets)
  mdofs_data  = Vector{Int32}(undef, n_offsets)
  coeffs_data = Vector{CT}(undef, n_offsets)
  ptrs = zeros(eltype(sDOF_to_mdofs.ptrs), n_slaves + 1)
  ptrs[1] = 1
  kf = 0
  for s in 1:n_slaves
    b = sDOF_to_offsets[s]
    ptrs[s+1] = ptrs[s]
    if !iszero(b)
      kf += 1
      dmdof_to_offsets[kf] = OT(b)
      mdofs_data[kf]  = Int32(-(n_dmdofs + kf))
      coeffs_data[kf] = one(CT)
      ptrs[s+1] += 1
    end
  end

  return (
    vcat(mDOF_to_dof, zeros(Int32, n_offsets)),
    append_tables_locally(sDOF_to_mdofs,  Table(mdofs_data,  ptrs)),
    append_tables_locally(sDOF_to_coeffs, Table(coeffs_data, ptrs)),
    dmdof_to_offsets
  )
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
  n_mdofs = length(mDOF_to_dof)
  n_fdofs = num_free_dofs(space)
  n_dofs  = n_fdofs + num_dirichlet_dofs(space)
  DOF_to_msDOF = Vector{Int}(undef,n_dofs)

  for (mDOF,dof) in enumerate(mDOF_to_dof)
    iszero(dof) && continue
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_to_msDOF[DOF] = mDOF
  end

  for (sDOF,dof) in enumerate(sDOF_to_dof)
    DOF = _dof_to_DOF(dof,n_fdofs)
    DOF_to_msDOF[DOF] = sDOF + n_mdofs
  end

  return DOF_to_msDOF
end

# ---------------------------------------------------------------------------
# Table-level merge and chain-resolution for constraint tables
# ---------------------------------------------------------------------------

"""
    merge_slave_constraint_tables(
        space, s1_dof, s1_dofs, s1_coeffs, s2_dof, s2_dofs, s2_coeffs,
        s1_offsets, s2_offsets; on_conflict)
        → (sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets)

Merge two Constructor-2 constraint table triples defined on the same `space`.
Each triple lists only slave DOFs: `s_dof[i]` is the signed dof of slave `i`,
`s_dofs[i]` its master signed dofs, `s_coeffs[i]` the corresponding coefficients.
`s1_offsets` / `s2_offsets` are optional per-slave inhomogeneity vectors (default zero).

A conflict arises when the same slave dof appears in both sets:
  - `on_conflict=nothing` (default): raise an error.
  - otherwise: `on_conflict(dof, dofs1, coeffs1, offset1, dofs2, coeffs2, offset2) → (dofs, coeffs, offset)`
    is called and its return value replaces the stored row.
"""
function merge_slave_constraint_tables(
  space::SingleFieldFESpace,
  s1_dof, s1_dofs::Table, s1_coeffs::Table{T},
  s2_dof, s2_dofs::Table, s2_coeffs::Table{T},
  s1_offsets::AbstractVector = zeros(T, length(s1_dof)),
  s2_offsets::AbstractVector = zeros(T, length(s2_dof));
  on_conflict=nothing
) where T
  n_fdofs = num_free_dofs(space)
  n_dofs  = n_fdofs + num_dirichlet_dofs(space)
  n_s1    = length(s1_dof)
  n_s2    = length(s2_dof)

  s1_DOF    = [_dof_to_DOF(dof, n_fdofs) for dof in s1_dof]
  DOF_to_s1 = find_inverse_index_map(s1_DOF, n_dofs)

  # conflicted_s1[i1] = i2 (s2 index of the conflicting entry) or 0 (no conflict).
  conflicted_s1  = zeros(Int, n_s1)
  new_s2_indices = Int[]
  sizehint!(new_s2_indices, n_s2)

  for i2 in 1:n_s2
    i1 = DOF_to_s1[_dof_to_DOF(s2_dof[i2], n_fdofs)]
    if iszero(i1)
      push!(new_s2_indices, i2)
    elseif isnothing(on_conflict)
      error("Constraint conflict on dof $(s2_dof[i2]):\n" *
            "  existing: $(collect(dataview(s1_dofs, i1)))\n" *
            "  incoming: $(collect(dataview(s2_dofs, i2)))")
    else
      conflicted_s1[i1] = i2
    end
  end

  # Fast path: no conflicts.
  if all(iszero, conflicted_s1)
    return vcat(s1_dof, s2_dof),
           append_tables_globally(s1_dofs,   s2_dofs),
           append_tables_globally(s1_coeffs, s2_coeffs),
           vcat(s1_offsets, s2_offsets)
  end

  n_new_s2 = length(new_s2_indices)
  n_out    = n_s1 + n_new_s2

  # Pessimistic data size: conflicts get len(s1_row) + len(s2_row) as upper bound.
  n_data = 0
  for i in 1:n_s1
    i2 = conflicted_s1[i]
    n_data += s1_dofs.ptrs[i+1] - s1_dofs.ptrs[i]
    !iszero(i2) && (n_data += s2_dofs.ptrs[i2+1] - s2_dofs.ptrs[i2])
  end
  for i2 in new_s2_indices
    n_data += s2_dofs.ptrs[i2+1] - s2_dofs.ptrs[i2]
  end

  out_dof     = Vector{Int}(undef, n_out)
  out_offsets = Vector{T}(undef, n_out)
  out_data_d  = Vector{Int}(undef, n_data)
  out_data_c  = Vector{T}(undef, n_data)
  ptrs        = zeros(Int32, n_out + 1)
  copyto!(out_dof, 1, s1_dof, 1, n_s1)

  k = 0
  for i1 in 1:n_s1
    i2 = conflicted_s1[i1]
    if iszero(i2)
      rin = datarange(s1_dofs, i1)
      nl  = length(rin)
      out_data_d[k+1:k+nl] = view(s1_dofs.data,   rin)
      out_data_c[k+1:k+nl] = view(s1_coeffs.data, rin)
      out_offsets[i1] = s1_offsets[i1]
    else
      r1 = datarange(s1_dofs, i1);  r2 = datarange(s2_dofs, i2)
      od, oc, oo = on_conflict(s1_dof[i1],
        view(s1_dofs.data, r1), view(s1_coeffs.data, r1), s1_offsets[i1],
        view(s2_dofs.data, r2), view(s2_coeffs.data, r2), s2_offsets[i2])
      nl = length(od)
      out_data_d[k+1:k+nl] = od
      out_data_c[k+1:k+nl] = oc
      out_offsets[i1] = oo
    end
    ptrs[i1+1] = nl
    k += nl
  end
  for (j, i2) in enumerate(new_s2_indices)
    ii2 = n_s1 + j
    rin = datarange(s2_dofs, i2)
    nl  = length(rin)
    out_dof[ii2]          = s2_dof[i2]
    out_offsets[ii2]      = s2_offsets[i2]
    out_data_d[k+1:k+nl]  = view(s2_dofs.data,   rin)
    out_data_c[k+1:k+nl]  = view(s2_coeffs.data, rin)
    ptrs[ii2+1] = nl
    k += nl
  end

  Arrays.length_to_ptrs!(ptrs)
  resize!(out_data_d, k)
  resize!(out_data_c, k)
  return out_dof, Table(out_data_d, ptrs), Table(out_data_c, ptrs), out_offsets
end

"""
    compute_topological_ordering(dag::Table{Int}, in_degree::Vector{Int}; keys=1:length(dag))
        → (order::Vector{Int}, order_rev::Vector{Int}, has_chains::Bool)

Core of the topological sort shared by `close_slave_constraint_tables` and
`ConstraintHandler`'s `close!`.

- `dag` is the "dependents" graph in CSR (`Table`) form: `dag[t]` contains the 
nodes that depend on node `t`. 
In the context of constraint resolution, an edge `t → s` means
that slave `s` depends on slave `t`.
- `in_degree[s]` is the number of (as yet unresolved) dependencies of `s`; it is mutated in place.

Kahn's algorithm: entities with `in_degree == 0` are processed in ascending
order of `keys[s]`, via a priority queue. This makes the resulting order
deterministic and independent of array/insertion order — `keys` defaults to
`1:n` (slave-array order) but can be set to e.g. global DoF ids.

`order[i]` is the slave (in `1:n`) at topological position `i`; `order_rev`
is its inverse (`order_rev[order[i]] == i`). `has_chains` is `false` iff `dag`
has no edges at all, i.e. no slave's master is itself a slave — in that case
`order`/`order_rev` carry no useful information beyond `keys`-sorting, and
callers can fast-path (resolution is a no-op).

Errors, listing the offending `keys`, if a cycle remains (i.e. not all `n`
entities could be ordered).
"""
function compute_topological_ordering(
  dag::Table{Int}, in_degree::Vector{Int};
  keys = 1:length(dag)
)
  n = length(dag)
  has_chains = !isempty(dag.data)

  order     = Vector{Int}(undef, n)
  order_rev = Vector{Int}(undef, n)

  pq = PriorityQueue{Int,eltype(keys)}()
  for s in 1:n
    iszero(in_degree[s]) && enqueue!(pq, s, keys[s])
  end

  head = 0
  while !isempty(pq)
    s = dequeue!(pq)
    head += 1
    order[head]  = s
    order_rev[s] = head
    for t in dataview(dag, s)
      in_degree[t] -= 1
      iszero(in_degree[t]) && enqueue!(pq, t, keys[t])
    end
  end

  if head < n
    cycle = sort!([keys[s] for s in 1:n if in_degree[s] > 0])
    error("Constraint cycle detected among DoFs: $cycle\n" *
          "Cycles are not resolvable by substitution. Check your constraint sources.")
  end

  return order, order_rev, has_chains
end

"""
    compute_topological_ordering(args...; kwargs...)
        → (order, order_rev, has_chains)

Ingests either a Constructor-2 constraint table
(`sDOF_to_dofs, n_fdofs, dof_to_sDOF`) or a `ConstraintHandler`, via
[`compute_constraint_dag`](@ref), and topologically sorts the resulting
slave→slave dependency graph. See the core method for details.
"""
function compute_topological_ordering(args...; kwargs...)
  dag, in_degree = compute_constraint_dag(args...)
  return compute_topological_ordering(dag, in_degree; kwargs...)
end

"""
    compute_constraint_dag(sDOF_to_dofs, n_fdofs, dof_to_sDOF)
        → (dag::Table{Int}, in_degree::Vector{Int})

Build the slave→slave dependency graph for a Constructor-2 constraint table:
slave `s` depends on slave `t` if `t`'s DoF appears as a master of `s`.
"""
function compute_constraint_dag(sDOF_to_dofs::Table, n_fdofs::Int, dof_to_sDOF::AbstractVector)
  n_slaves = length(sDOF_to_dofs)

  dep_ptrs  = zeros(Int32, n_slaves + 1)
  in_degree = zeros(Int, n_slaves)
  for s in 1:n_slaves
    for d in dataview(sDOF_to_dofs, s)
      t = dof_to_sDOF[_dof_to_DOF(d, n_fdofs)]
      if t != 0
        dep_ptrs[t+1] += Int32(1)
        in_degree[s]  += 1
      end
    end
  end

  Arrays.length_to_ptrs!(dep_ptrs)
  dep_data = Vector{Int}(undef, dep_ptrs[end] - 1)
  fill_pos = copy(dep_ptrs)
  for s in 1:n_slaves
    for d in dataview(sDOF_to_dofs, s)
      t = dof_to_sDOF[_dof_to_DOF(d, n_fdofs)]
      if t != 0
        dep_data[fill_pos[t]] = s
        fill_pos[t] += 1
      end
    end
  end

  dag = Table(dep_data, dep_ptrs)
  return dag, in_degree
end

"""
    close_slave_constraint_tables(
        space, sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets; tol)
        → (sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets)

Chain-resolve a Constructor-2 constraint table triple: substitute any slave master
with its own constraint row until all masters are free. Errors on cycles.
Merges duplicate masters and strips near-zero coefficients (threshold `tol`).
The slave set is invariant under resolution; output rows are in topological order.

`sDOF_to_offsets` (optional, default zero) carries inhomogeneities: if slave i
depends on slave j the offset accumulates as `offset_i += coeff_ij * offset_j`.
"""
function close_slave_constraint_tables(
  space::SingleFieldFESpace,
  sDOF_to_dof, sDOF_to_dofs::Table, sDOF_to_coeffs::Table{T},
  sDOF_to_offsets::AbstractVector = zeros(T, length(sDOF_to_dof));
  tol::T = 1000*eps(T), keys = sDOF_to_dof
) where T
  n_fdofs  = num_free_dofs(space)
  n_dofs   = n_fdofs + num_dirichlet_dofs(space)
  n_slaves = length(sDOF_to_dof)

  s_DOF       = [_dof_to_DOF(dof, n_fdofs) for dof in sDOF_to_dof]
  dof_to_sDOF = find_inverse_index_map(s_DOF, n_dofs)

  order, order_rev, has_chains = compute_topological_ordering(sDOF_to_dofs, n_fdofs, dof_to_sDOF; keys)
  has_chains || return sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, copy(sDOF_to_offsets)

  # Exact pre-dedup/strip output size: single forward pass in topo order.
  n_data = 0
  ub     = zeros(Int, n_slaves)
  for i in 1:n_slaves
    s = order[i]
    for d in dataview(sDOF_to_dofs, s)
      t = dof_to_sDOF[_dof_to_DOF(d, n_fdofs)]
      ub[i] += (t != 0) ? ub[order_rev[t]] : 1
    end
    n_data += ub[i]
  end

  out_offsets = copy(sDOF_to_offsets)
  out_data_d  = Vector{Int}(undef, n_data)
  out_data_c  = Vector{T}(undef, n_data)
  ptrs        = zeros(Int32, n_slaves + 1)
  ptrs[1]     = Int32(1)
  acc         = Dict{Int, T}()
  k           = 0

  for i in 1:n_slaves
    s  = order[i]
    rd = dataview(sDOF_to_dofs, s)
    rc = dataview(sDOF_to_coeffs, s)

    any_sub = any(dof_to_sDOF[_dof_to_DOF(d, n_fdofs)] != 0 for d in rd)

    if any_sub
      empty!(acc)
      for kk in eachindex(rd)
        d = rd[kk];  c = rc[kk]
        t = dof_to_sDOF[_dof_to_DOF(d, n_fdofs)]
        if t != 0
          j = order_rev[t]
          for m in ptrs[j]:ptrs[j+1]-1
            dd = out_data_d[m]
            acc[dd] = get(acc, dd, zero(T)) + c * out_data_c[m]
          end
          out_offsets[s] += c * out_offsets[t]
        else
          acc[d] = get(acc, d, zero(T)) + c
        end
      end
    end

    for (d, c) in (any_sub ? acc : zip(rd, rc))
      abs(c) > tol || continue
      k += 1;  out_data_d[k] = d;  out_data_c[k] = c
    end
    ptrs[i+1] = Int32(k + 1)
  end

  resize!(out_data_d, k)
  resize!(out_data_c, k)
  out_dof = [sDOF_to_dof[order[i]]    for i in 1:n_slaves]
  out_off = [out_offsets[order[i]] for i in 1:n_slaves]
  return out_dof, Table(out_data_d, ptrs), Table(out_data_c, ptrs), out_off
end
