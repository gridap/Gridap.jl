######################
# ODEOpFromFEOpCache #
######################
"""
    struct ODEOpFromFEOpCache <: GridapType

Structure that stores the `TransientFESpace` and cache of a
`TransientFEOperator`, as well as the jacobian matrices and residual if they
are constant.
"""
mutable struct ODEOpFromFEOpCache <: GridapType
  Us
  Uts
  feopcache
  jac_const_mats
  jacvec
  res_const_vec
end

#################
# ODEOpFromFEOp #
#################
"""
    struct ODEOpFromFEOp <: ODEOperator end

Wrapper that transforms a `TransientFEOperator` into an `ODEOperator`, i.e.
takes `residual(t, uh, ∂t[uh], ..., ∂t^N[uh], vh)` and returns
`residual(t, us)`, where `us[k] = ∂t^k[us]` and `uf` represents the free values
of `uh`.
"""
struct ODEOpFromFEOp{C} <: ODEOperator{C}
  feop::TransientFEOperator{C}
end

# ODEOperator interface
Polynomials.get_order(odeop::ODEOpFromFEOp) = get_order(odeop.feop)

function allocate_odeopcache(
  odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}}
)
  # Allocate FE spaces for derivatives
  order = get_order(odeop)
  Ut = get_trial(odeop.feop)
  U = allocate_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for k in 1:order
    Uts = (Uts..., ∂t(Uts[k]))
    Us = (Us..., allocate_space(Uts[k+1]))
  end

  # Allocate the cache of the FE operator
  feopcache = allocate_feopcache(odeop.feop)

  # Create uh from us for assembly
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  uh = _make_uh_from_us(odeop, us, Us)

  # Store the jacobians that are constant
  jac_const_mats = ()
  use_jacvec = false
  for k in 0:order
    jac_const_mat = nothing
    if is_jacobian_constant(odeop, k)
      matdata = _matdata_jacobian(odeop.feop, t, uh, k, 1)
      jac_const_mat = assemble_matrix(get_assembler(odeop.feop), matdata)
      use_jacvec = true
    end
    jac_const_mats = (jac_const_mats..., jac_const_mat)
  end

  # Allocate a vector to compute the jacobian vector product if any jacobian
  # has been stored
  jacvec = nothing
  if use_jacvec
    res = get_res(odeop.feop)
    vecdata = collect_cell_vector(V, res(t, uh, v))
    jacvec = allocate_vector(get_assembler(odeop.feop), vecdata)
  end

  # Store the forcing term if it is constant
  res_const_vec = nothing
  if is_residual_constant(odeop)
    res = get_res(odeop.feop)
    vecdata = collect_cell_vector(V, res(t, uh, v))
    res_const_vec = assemble_vector(get_assembler(odeop.feop), vecdata)
  end

  ODEOpFromFEOpCache(Us, Uts, feopcache, jac_const_mats, jacvec, res_const_vec)
end

function update_odeopcache!(odeopcache, odeop::ODEOpFromFEOp, t::Real)
  Us = ()
  for k in 0:get_order(odeop)
    Us = (Us..., evaluate!(odeopcache.Us[k+1], odeopcache.Uts[k+1], t))
  end
  odeopcache.Us = Us
  odeopcache.feopcache = update_feopcache!(odeopcache.feopcache, odeop.feop, t)
  odeopcache
end

function Algebra.allocate_residual(
  odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  res = get_res(odeop.feop)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  allocate_vector(get_assembler(odeop.feop), vecdata)
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache;
  filter::Tuple{Vararg{Bool}}=ntuple(_ -> true, get_order(odeop) + 2)
)
  if any(filter)
    uh = _make_uh_from_us(odeop, us, odeopcache.Us)
    V = get_test(odeop.feop)
    v = get_fe_basis(V)
    res = get_res(odeop.feop)
    vecdata = collect_cell_vector(V, res(t, uh, v))
    assemble_vector!(r, get_assembler(odeop.feop), vecdata)
  else
    fill!(r, zero(eltype(r)))
  end

  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromFEOp{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache;
  filter::Tuple{Vararg{Bool}}=ntuple(_ -> true, get_order(odeop) + 2)
)
  order = get_order(odeop)
  res_needed = any(i -> filter[i], 1:order+1)
  res_const = !isnothing(odeopcache.res_const_vec)

  mass_needed = filter[order+2]
  mass_const = !isnothing(odeopcache.jac_const_mats[end])

  # Compute uh from us if assembly is needed
  if (res_needed && !res_const) || (mass_needed && !mass_const)
    V = get_test(odeop.feop)
    v = get_fe_basis(V)
    uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  end

  fill!(r, zero(eltype(r)))

  # Add residual (without mass) if needed
  if res_needed
    if res_const
      axpy!(1, odeopcache.res_const_vec, r)
    else
      res = get_res(odeop.feop)
      vecdata = collect_cell_vector(V, res(t, uh, v))
      assemble_vector_add!(r, get_assembler(odeop.feop), vecdata)
    end
  end

  # Add mass if needed
  if mass_needed
    if mass_const
      mul!(odeopcache.jacvec, odeopcache.jac_const_mats[end], us[end])
      axpy!(1, odeopcache.jacvec, r)
    else
      # Build the time derivative
      ∂tNuh = uh
      for _ in 1:order
        ∂tNuh = ∂t(∂tNuh)
      end
      mass = get_mass(odeop.feop)
      if odeop isa AbstractSemilinearODE
        vecdata = collect_cell_vector(V, mass(t, ∂tNuh, v))
      else
        vecdata = collect_cell_vector(V, mass(t, uh, ∂tNuh, v))
      end
      assemble_vector_add!(r, get_assembler(odeop.feop), vecdata)
    end
  end

  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache;
  filter::Tuple{Vararg{Bool}}=ntuple(_ -> true, get_order(odeop) + 2)
)
  order = get_order(odeop)
  res_needed = filter[1]
  res_const = !isnothing(odeopcache.res_const_vec)

  form_needed_not_const = any(
    k -> filter[k+2] && isnothing(odeopcache.jac_const_mats[k+1]),
    0:order
  )

  # Compute uh from us if assembly is needed
  if (res_needed && !res_const) || form_needed_not_const
    V = get_test(odeop.feop)
    v = get_fe_basis(V)
    uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  end

  fill!(r, zero(eltype(r)))

  # Add the forcing term if needed
  if res_needed
    if res_const
      axpy!(1, odeopcache.res_const_vec, r)
    else
      res = get_res(odeop.feop)
      vecdata = collect_cell_vector(V, res(t, uh, v))
      assemble_vector_add!(r, get_assembler(odeop.feop), vecdata)
    end
  end

  # Add the forms if needed
  ∂tkuh = uh
  for k in 0:order
    idx_form = k + 1
    idx_filter = k + 2
    if filter[idx_filter]
      jac_mat = odeopcache.jac_const_mats[idx_form]
      form_const = !isnothing(jac_mat)
      if form_const
        mul!(odeopcache.jacvec, jac_mat, us[idx_form])
        axpy!(1, odeopcache.jacvec, r)
      else
        form = get_forms(odeop.feop)[idx_form]
        println(get_forms(odeop.feop))
        vecdata = collect_cell_vector(V, form(t, ∂tkuh, v))
        assemble_vector_add!(r, get_assembler(odeop.feop), vecdata)
      end
    end
    if k < order
      ∂tkuh = ∂t(∂tkuh)
    end
  end

  r
end

function Algebra.allocate_jacobian(
  odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  _matdata = ()
  for k in 0:get_order(odeop.feop)
    _matdata = (_matdata..., _matdata_jacobian(odeop.feop, t, uh, k, 0))
  end
  matdata = _vcat_matdata(_matdata)
  allocate_matrix(get_assembler(odeop.feop), matdata)
end

function Algebra.jacobian!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  # Retrieve stored jacobian or compute and assemble the domain contribution
  if !isnothing(odeopcache.jac_const_mats[k+1])
    axpy_sparse!(γ, odeopcache.jac_const_mats[k+1], J)
  else
    uh = _make_uh_from_us(odeop, us, odeopcache.Us)
    matdata = _matdata_jacobian(odeop.feop, t, uh, k, γ)
    assemble_matrix_add!(J, get_assembler(odeop.feop), matdata)
  end

  J
end

function jacobians!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  odeopcache
)
  need_assembly = any(isnothing, odeopcache.jac_const_mats)
  if need_assembly
    uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  end

  # Loop through jacobians and check whether they are already stored,
  # otherwise compute and add the domain contributions
  _matdata = ()
  for k in 0:get_order(odeop)
    γ = γs[k+1]
    if !iszero(γ)
      if !isnothing(odeopcache.jac_const_mats[k+1])
        axpy_sparse!(γ, odeopcache.jac_const_mats[k+1], J)
      else
        _matdata = (_matdata..., _matdata_jacobian(odeop.feop, t, uh, k, γ))
      end
    end
  end

  # Assemble the domain contributions if needed
  if need_assembly
    matdata = _vcat_matdata(_matdata)
    assemble_matrix_add!(J, get_assembler(odeop.feop), matdata)
  end

  J
end

function is_jacobian_constant(odeop::ODEOpFromFEOp, k::Integer)
  is_jacobian_constant(odeop.feop, k)
end

function is_residual_constant(odeop::ODEOpFromFEOp)
  is_residual_constant(odeop.feop)
end

#########
# Utils #
#########
function _matdata_jacobian(
  odeop::TransientFEOperatorTypes,
  t::Real, uh::TransientCellField,
  k::Integer, γ::Real
)
  Ut = evaluate(get_trial(odeop), nothing)
  V = get_test(odeop)
  du = get_trial_fe_basis(Ut)
  v = get_fe_basis(V)
  jac = get_jacs(odeop)[k+1]
  collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
end

function _vcat_matdata(_matdata)
  term_to_cellmat_j = ()
  term_to_cellidsrows_j = ()
  term_to_cellidscols_j = ()
  for j in 1:length(_matdata)
    term_to_cellmat_j = (term_to_cellmat_j..., _matdata[j][1])
    term_to_cellidsrows_j = (term_to_cellidsrows_j..., _matdata[j][2])
    term_to_cellidscols_j = (term_to_cellidscols_j..., _matdata[j][3])
  end
  term_to_cellmat = vcat(term_to_cellmat_j...)
  term_to_cellidsrows = vcat(term_to_cellidsrows_j...)
  term_to_cellidscols = vcat(term_to_cellidscols_j...)
  (term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
end

function _make_uh_from_us(odeop, us, Us)
  u = EvaluationFunction(Us[1], us[1])
  dus = ()
  for k in 1:get_order(odeop)
    dus = (dus..., EvaluationFunction(Us[k+1], us[k+1]))
  end
  TransientCellField(u, dus)
end
