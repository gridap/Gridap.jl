
struct ZeroMeanFESpace{CA,S,D<:AbstractVector{Float64}} <: SingleFieldFESpace
  space::FESpaceWithConstantFixed{CA,S}
  int_φi::Vector{Float64}
  vol::Float64
  h_free::D
  h_diri::Float64
end

"""
    ZeroMeanFESpace(U::SingleFieldFESpace, dΩ::Measure)

Subspace of `U` of the functions with zero mean, that is all `u ∈ U`
such that

    ∫(u)dΩ = 0

One DOF of `U` (whose shape-function has non-zera mean) is fixed to force the
constraint.

# Extended help

In order to be robust to non-partition of unity dof basis of `U`, the space is
"decomposed" as

    u = Σ_i u_i φ_i = (Σ_i≠k f_i φ_i) + f_k φ_k + c 𝔥

where `k` is the DOF of `U` that the ZeroMeanFESpace fixes, `u_i = f_i +
c·h_i`, `φ_i` are `U`'s shape functions, and `𝔥 ∈ U` is the Riesz representative
of the mean functional,

    ⟨𝔥,q⟩ = ∫(𝔥*q)dΩ = ∫q  for all q.

`𝔥`'s dof values, `h_dofs`, is the vector of ones for a basis that is a
partition of unity (like Lagrangian).

So, `u` is imposed to have zero mean by fixing the constant `c` to

c = - ( (Σ_i≠k f_i ∫φ_i) + d ∫φ_k) / "vol".

Evaluating `u` still requires `h_dofs`.
"""
function ZeroMeanFESpace(space::SingleFieldFESpace, dΩ::Measure)
  int_φi = assemble_vector(v->∫(v)*dΩ, space)
  T = eltype(int_φi)
  @notimplementedif !isone(num_indep_components(T))
  uh = interpolate(one(T), space)
  h_dofs = get_free_dof_values(uh)
  h_dofs = all(isone,h_dofs) ? Ones(length(h_dofs)) : h_dofs
  vol = h_dofs isa Ones ? sum(int_φi) : dot(int_φi, h_dofs)

  vol > zero(vol) || throw(ArgumentError("""\n
  A `ZeroMeanFESpace` of volume $vol has no function of unit mean value to
  shift its functions along, so no function of it can be given zero mean value.
  """))
  @check _represents_mean_functional(space, dΩ, int_φi, uh) """\n
  `space` does not contain a riesz representative for the mean functional.
  Either it does not contain the constants and you should impose the zero mean
  value with a Lagrange multiplier (using a mixed field formulation with a
  `ConstantFESpace` unknown), or the space's DOF basis is broken.
  """

  if h_dofs isa Ones
    dof_to_fix = length(h_dofs)
    h_free = Ones(length(h_dofs)-1)
  else
    # The h_dofs of the dof to fix must not be zero
    dof_to_fix = last(findmax(abs, h_dofs))
    h_free = deleteat!(copy(h_dofs), dof_to_fix)
  end

  _space = FESpaceWithConstantFixed(space, true, dof_to_fix)
  ZeroMeanFESpace(_space, int_φi, vol, h_free, h_dofs[dof_to_fix])
end

# Checks whether uh is the Riesz representative of the mean functional for the
# φ -> ∫(uh * φ) inner product, that is that q -> ∫( q(𝔥-1) ) is zero.
function _represents_mean_functional(space, dΩ, int_φi, uh)
  one_ = one(eltype(int_φi))
  norm(assemble_vector(q -> ∫( *(q, (uh-one_)) )*dΩ,space)) ≤ 1.0e-6*norm(int_φi)
end

function TrialFESpace(f::ZeroMeanFESpace)
  U = TrialFESpace(f.space)
  ZeroMeanFESpace(U, f.int_φi, f.vol, f.h_free, f.h_diri)
end

# function FEFunction(f::ZeroMeanFESpace,free_values)
#   msg = """
#     This function should never be called for ZeroMeanFESpace. Please use
#     `FEFunction(f::ZeroMeanFESpace,free_values,dirichlet_values)` instead.
#     Reason:
#       Without the fixed value (i.e the dirichlet_values), we cannot correctly
#       interpolate the free dofs within the space.
#   """
#   @unreachable msg
# end
#
function FEFunction(
  f::ZeroMeanFESpace,
  free_values::AbstractVector,
  dirichlet_values::AbstractVector # only 1, the average fixing dof
)
  c = _compute_new_fixedval(
    free_values,
    dirichlet_values,
    f.int_φi,
    f.vol,
    f.space.dof_to_fix
  )
  fv = lazy_map(+, free_values, c .* f.h_free)
  dv = dirichlet_values .+ c*f.h_diri
  FEFunction(f.space, fv, dv)
end

function EvaluationFunction(f::ZeroMeanFESpace,free_values)
  FEFunction(f.space,free_values)
end

function _compute_new_fixedval(fv, dv, int_φi, vol, fixed_dof)
  @assert length(fv) + 1 == length(int_φi)
  @assert length(dv) == 1
  c = zero(eltype(int_φi))
  for i=1:fixed_dof-1
    c += fv[i]*int_φi[i]
  end
  c += first(dv)*int_φi[fixed_dof]
  for i=fixed_dof+1:length(int_φi)
    c += fv[i-1]*int_φi[i]
  end
  c = -c/vol
  return c
end

# This is required, otherwise we end up calling `FEFunction` with a fixed value of zero,
# which does not properly interpolate the function provided.
# With this change, we are interpolating in the unconstrained space and then
# substracting the mean.
function interpolate!(object,free_values,fs::ZeroMeanFESpace)
  dirichlet_values = zero_dirichlet_values(fs)
  interpolate_everywhere!(object,free_values,dirichlet_values,fs)
end

# Delegated functions

get_triangulation(f::ZeroMeanFESpace) = get_triangulation(f.space)

ConstraintStyle(::Type{<:ZeroMeanFESpace{CA,S}}) where {CA,S} = ConstraintStyle(S)

CellField(t::ZeroMeanFESpace,cell_vals) = CellField(t.space,cell_vals)

get_cell_isconstrained(f::ZeroMeanFESpace) = get_cell_isconstrained(f.space)

get_cell_constraints(f::ZeroMeanFESpace) = get_cell_constraints(f.space)

get_dirichlet_dof_values(f::ZeroMeanFESpace) = get_dirichlet_dof_values(f.space)

get_fe_basis(f::ZeroMeanFESpace) = get_fe_basis(f.space)

get_trial_fe_basis(f::ZeroMeanFESpace) = get_trial_fe_basis(f.space)

get_fe_dof_basis(f::ZeroMeanFESpace) = get_fe_dof_basis(f.space)

get_free_dof_ids(f::ZeroMeanFESpace) = get_free_dof_ids(f.space)

zero_free_values(f::ZeroMeanFESpace) = zero_free_values(f.space)

get_vector_type(f::ZeroMeanFESpace) = get_vector_type(f.space)

get_dof_value_type(f::ZeroMeanFESpace) = get_dof_value_type(f.space)

get_cell_dof_ids(f::ZeroMeanFESpace) = get_cell_dof_ids(f.space)

get_dirichlet_dof_ids(f::ZeroMeanFESpace) = get_dirichlet_dof_ids(f.space)

zero_dirichlet_values(f::ZeroMeanFESpace) = zero_dirichlet_values(f.space)

num_dirichlet_tags(f::ZeroMeanFESpace) = num_dirichlet_tags(f.space)

get_dirichlet_dof_tag(f::ZeroMeanFESpace) = get_dirichlet_dof_tag(f.space)

scatter_free_and_dirichlet_values(f::ZeroMeanFESpace,fv,dv) = scatter_free_and_dirichlet_values(f.space,fv,dv)

gather_free_and_dirichlet_values(f::ZeroMeanFESpace,cv) = gather_free_and_dirichlet_values(f.space,cv)

gather_free_and_dirichlet_values!(fv,dv,f::ZeroMeanFESpace,cv) = gather_free_and_dirichlet_values!(fv,dv,f.space,cv)

gather_dirichlet_values(f::ZeroMeanFESpace,cv) = gather_dirichlet_values(f.space,cv)

gather_dirichlet_values!(dv,f::ZeroMeanFESpace,cv) = gather_dirichlet_values!(dv,f.space,cv)

gather_free_values(f::ZeroMeanFESpace,cv) = gather_free_values(f.space,cv)

gather_free_values!(fv,f::ZeroMeanFESpace,cv) = gather_free_values!(fv,f.space,cv)
