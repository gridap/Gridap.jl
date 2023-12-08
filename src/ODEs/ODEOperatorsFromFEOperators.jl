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
  const_forms
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

function get_forms(odeop::ODEOpFromFEOp)
  get_forms(odeop.feop)
end

function is_form_constant(odeop::ODEOpFromFEOp, k::Integer)
  is_form_constant(odeop.feop, k)
end

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
  feopcache = allocate_feopcache(odeop.feop, t, us)

  # Variables for assembly
  uh = _make_uh_from_us(odeop, us, Us)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  Ut = evaluate(get_trial(odeop.feop), nothing)
  du = get_trial_fe_basis(Ut)
  assembler = get_assembler(odeop.feop)

  # Store the forms that are constant
  const_forms = ()
  forms = get_forms(odeop.feop)
  jacs = get_jacs(odeop.feop)

  # Little workaround here since when the `ODEOperator` is quasilinear or
  # semilinear but not linear, it has only one form but `order+1` jacobians
  # so we need to be careful with indexing
  odeoptype = ODEOperatorType(odeop)
  if odeoptype <: AbstractLinearODE
    for k in 0:length(forms)-1
      const_form = nothing
      if is_form_constant(odeop, k)
        jac = jacs[k+1]
        matdata = collect_cell_matrix(Ut, V, jac(t, uh, du, v))
        const_form = assemble_matrix(assembler, matdata)
      end
      const_forms = (const_forms..., const_form)
    end
  elseif odeoptype <: AbstractQuasilinearODE
    const_form = nothing
    k = order
    if is_form_constant(odeop, k)
      jac = jacs[k+1]
      matdata = collect_cell_matrix(Ut, V, jac(t, uh, du, v))
      const_form = assemble_matrix(assembler, matdata)
    end
    const_forms = (const_forms..., const_form)
  end

  ODEOpFromFEOpCache(Us, Uts, feopcache, const_forms)
end

function update_odeopcache!(odeopcache, odeop::ODEOpFromFEOp, t::Real)
  Us = ()
  for k in 0:get_order(odeop)
    Us = (Us..., evaluate!(odeopcache.Us[k+1], odeopcache.Uts[k+1], t))
  end
  odeopcache.Us = Us

  feopcache, feop = odeopcache.feopcache, odeop.feop
  odeopcache.feopcache = update_feopcache!(feopcache, feop, t)

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
  assembler = get_assembler(odeop.feop)

  res = get_res(odeop.feop)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  allocate_vector(assembler, vecdata)
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  res = get_res(odeop.feop)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  assemble_vector!(r, assembler, vecdata)
  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromFEOp{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  fill!(r, zero(eltype(r)))
  # Mass
  order = get_order(odeop)
  mass = get_forms(odeop.feop)[1]
  ∂tNuh = ∂t(uh, Val(order))
  vecdata = collect_cell_vector(V, mass(t, uh, ∂tNuh, v))
  assemble_vector_add!(r, assembler, vecdata)

  # Residual
  res = get_res(odeop.feop)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromFEOp{<:AbstractSemilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  fill!(r, zero(eltype(r)))
  # Mass
  order = get_order(odeop)
  mass = get_forms(odeop.feop)[1]
  ∂tNuh = ∂t(uh, Val(order))
  vecdata = collect_cell_vector(V, mass(t, ∂tNuh, v))
  assemble_vector_add!(r, assembler, vecdata)

  # Residual
  res = get_res(odeop.feop)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  fill!(r, zero(eltype(r)))
  # Forms
  order = get_order(odeop)
  forms = get_forms(odeop.feop)
  ∂tkuh = uh
  for k in 0:order
    form = forms[k+1]
    vecdata = collect_cell_vector(V, form(t, ∂tkuh, v))
    assemble_vector_add!(r, assembler, vecdata)
    if k < order
      ∂tkuh = ∂t(∂tkuh)
    end
  end

  # Residual
  res = get_res(odeop.feop)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Algebra.allocate_jacobian(
  odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.feop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  all_matdata = ()
  jacs = get_jacs(odeop.feop)
  for k in 0:get_order(odeop.feop)
    jac = jacs[k+1]
    one_matdata = collect_cell_matrix(Ut, V, jac(t, uh, du, v))
    all_matdata = (all_matdata..., one_matdata)
  end
  matdata = _vcat_matdata(all_matdata)
  allocate_matrix(assembler, matdata)
end

function Algebra.jacobian!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.feop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  jac = get_jacs(odeop.feop)[k+1]
  matdata = collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
  assemble_matrix_add!(J, assembler, matdata)
  J
end

function Algebra.jacobian!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp{<:AbstractSemilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  order = get_order(odeop)
  # Special case for the mass matrix
  if k == order && is_form_constant(odeop, order)
    axpy_sparse!(γ, odeopcache.const_forms[1], J)
  else
    uh = _make_uh_from_us(odeop, us, odeopcache.Us)
    Ut = evaluate(get_trial(odeop.feop), nothing)
    du = get_trial_fe_basis(Ut)
    V = get_test(odeop.feop)
    v = get_fe_basis(V)
    assembler = get_assembler(odeop.feop)

    jac = get_jacs(odeop.feop)[k+1]
    matdata = collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
    assemble_matrix_add!(J, assembler, matdata)
  end
  J
end

function Algebra.jacobian!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  if is_form_constant(odeop, k)
    axpy_sparse!(γ, odeopcache.const_forms[k+1], J)
  else
    uh = _make_uh_from_us(odeop, us, odeopcache.Us)
    Ut = evaluate(get_trial(odeop.feop), nothing)
    du = get_trial_fe_basis(Ut)
    V = get_test(odeop.feop)
    v = get_fe_basis(V)
    assembler = get_assembler(odeop.feop)

    jac = get_jacs(odeop.feop)[k+1]
    matdata = collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
    assemble_matrix_add!(J, assembler, matdata)
  end
  J
end

function jacobians!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.feop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  all_matdata = ()
  jacs = get_jacs(odeop.feop)
  for k in 0:get_order(odeop)
    γ = γs[k+1]
    iszero(γ) && continue
    jac = jacs[k+1]
    one_matdata = collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
    all_matdata = (all_matdata..., one_matdata)
  end
  matdata = _vcat_matdata(all_matdata)
  assemble_matrix_add!(J, assembler, matdata)
  J
end

function jacobians!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.feop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  all_matdata = ()
  order = get_order(odeop)
  jacs = get_jacs(odeop.feop)
  for k in 0:order-1
    γ = γs[k+1]
    iszero(γ) && continue
    jac = jacs[k+1]
    one_matdata = collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
    all_matdata = (all_matdata..., one_matdata)
  end
  # Special case for the mass matrix
  k = order
  γ = γs[k+1]
  if !iszero(γ)
    if is_form_constant(odeop, k)
      axpy_sparse!(γ, odeopcache.const_forms[1], J)
    else
      jac = jacs[k+1]
      one_matdata = collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
      all_matdata = (all_matdata..., one_matdata)
    end
  end
  matdata = _vcat_matdata(all_matdata)
  assemble_matrix_add!(J, assembler, matdata)
  J
end

function jacobians!(
  J::AbstractMatrix, odeop::ODEOpFromFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.feop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.feop)

  all_matdata = ()
  jacs = get_jacs(odeop.feop)
  for k in 0:get_order(odeop)
    γ = γs[k+1]
    iszero(γ) && continue
    if is_form_constant(odeop, k)
      axpy_sparse!(γ, odeopcache.const_forms[k+1], J)
    else
      jac = jacs[k+1]
      one_matdata = collect_cell_matrix(Ut, V, γ * jac(t, uh, du, v))
      all_matdata = (all_matdata..., one_matdata)
    end
  end
  matdata = _vcat_matdata(all_matdata)
  assemble_matrix_add!(J, assembler, matdata)
  J
end

#########
# Utils #
#########
function _make_uh_from_us(odeop, us, Us)
  u = EvaluationFunction(Us[1], us[1])
  dus = ()
  for k in 1:get_order(odeop)
    dus = (dus..., EvaluationFunction(Us[k+1], us[k+1]))
  end
  TransientCellField(u, dus)
end

function _vcat_matdata(all_matdata)
  all_mats = ()
  all_rows = ()
  all_cols = ()
  for one_matdata in all_matdata
    one_mats, one_rows, one_cols = one_matdata
    all_mats = (all_mats..., one_mats)
    all_rows = (all_rows..., one_rows)
    all_cols = (all_cols..., one_cols)
  end
  flat_mats = vcat(all_mats...)
  flat_rows = vcat(all_rows...)
  flat_cols = vcat(all_cols...)
  (flat_mats, flat_rows, flat_cols)
end
