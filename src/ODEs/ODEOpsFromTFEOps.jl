#######################
# ODEOpFromTFEOpCache #
#######################
"""
    struct ODEOpFromTFEOpCache <: GridapType

Structure that stores the `TransientFESpace` and cache of a
`TransientFEOperator`, as well as the jacobian matrices and residual if they
are constant.
"""
mutable struct ODEOpFromTFEOpCache <: GridapType
  Us
  Uts
  tfeopcache
  const_forms
end

##################
# ODEOpFromTFEOp #
##################
"""
    struct ODEOpFromTFEOp <: ODEOperator end

Wrapper that transforms a `TransientFEOperator` into an `ODEOperator`, i.e.
takes `residual(t, uh, ∂t[uh], ..., ∂t^N[uh], vh)` and returns
`residual(t, us)`, where `us[k] = ∂t^k[us]` and `uf` represents the free values
of `uh`.
"""
struct ODEOpFromTFEOp{T} <: ODEOperator{T}
  tfeop::TransientFEOperator{T}

  function ODEOpFromTFEOp(tfeop::TransientFEOperator{T}) where {T}
    order = get_order(tfeop)
    if order == 0
      is_quasilinear = T <: AbstractQuasilinearODE
      is_linear = T <: AbstractLinearODE
      if is_quasilinear && !is_linear
        msg = """
        For an operator of order zero, the definitions of quasilinear,
        semilinear and linear coincide. Make sure that you have defined the
        transient FE operator as linear.
        """
        @unreachable msg
      else
        new{T}(tfeop)
      end
    else
      new{T}(tfeop)
    end
  end
end

# ODEOperator interface
function Polynomials.get_order(odeop::ODEOpFromTFEOp)
  get_order(odeop.tfeop)
end

function get_num_forms(odeop::ODEOpFromTFEOp)
  get_num_forms(odeop.tfeop)
end

function get_forms(odeop::ODEOpFromTFEOp)
  get_forms(odeop.tfeop)
end

function is_form_constant(odeop::ODEOpFromTFEOp, k::Integer)
  is_form_constant(odeop.tfeop, k)
end

function allocate_odeopcache(
  odeop::ODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}}
)
  # Allocate FE spaces for derivatives
  order = get_order(odeop)
  Ut = get_trial(odeop.tfeop)
  U = allocate_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for k in 1:order
    Uts = (Uts..., ∂t(Uts[k]))
    Us = (Us..., allocate_space(Uts[k+1]))
  end

  # Allocate the cache of the FE operator
  tfeopcache = allocate_tfeopcache(odeop.tfeop, t, us)

  # Variables for assembly
  uh = _make_uh_from_us(odeop, us, Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  assembler = get_assembler(odeop.tfeop)

  # Store the forms that are constant
  const_forms = ()
  num_forms = get_num_forms(odeop.tfeop)
  jacs = get_jacs(odeop.tfeop)

  # We want the stored jacobians to have the same sparsity as the full jacobian
  # (when all orders are considered), so we start by allocating it and we will assemble
  # the constant jacobians in a copy of the full jacobian
  # We need a little workaround here since when the `ODEOperator` is quasilinear or
  # semilinear but not linear, it has only one form but `order+1` jacobians.
  dc = DomainContribution()
  for k in 0:order
    jac = jacs[k+1]
    dc = dc + jac(t, uh, du, v)
  end
  matdata = collect_cell_matrix(Ut, V, dc)
  J_full = allocate_matrix(assembler, matdata)

  odeoptype = ODEOperatorType(odeop)
  if odeoptype <: AbstractLinearODE
    for k in 0:num_forms-1
      const_form = nothing
      if is_form_constant(odeop, k)
        jac = jacs[k+1]
        dc = jac(t, uh, du, v)
        matdata = collect_cell_matrix(Ut, V, dc)
        const_form = copy(J_full)
        fillstored!(const_form, zero(eltype(const_form)))
        assemble_matrix_add!(const_form, assembler, matdata)
      end
      const_forms = (const_forms..., const_form)
    end
  elseif odeoptype <: AbstractQuasilinearODE
    const_form = nothing
    k = order
    if is_form_constant(odeop, k)
      jac = jacs[k+1]
      dc = jac(t, uh, du, v)
      matdata = collect_cell_matrix(Ut, V, dc)
      const_form = copy(J_full)
      fillstored!(const_form, zero(eltype(const_form)))
      assemble_matrix_add!(const_form, assembler, matdata)
    end
    const_forms = (const_forms..., const_form)
  end

  ODEOpFromTFEOpCache(Us, Uts, tfeopcache, const_forms)
end

function update_odeopcache!(odeopcache, odeop::ODEOpFromTFEOp, t::Real)
  Us = ()
  for k in 0:get_order(odeop)
    Us = (Us..., evaluate!(odeopcache.Us[k+1], odeopcache.Uts[k+1], t))
  end
  odeopcache.Us = Us

  tfeopcache, tfeop = odeopcache.tfeopcache, odeop.tfeop
  odeopcache.tfeopcache = update_tfeopcache!(tfeopcache, tfeop, t)

  odeopcache
end

function Algebra.allocate_residual(
  odeop::ODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  res = get_res(odeop.tfeop)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  allocate_vector(assembler, vecdata)
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  res = get_res(odeop.tfeop)
  dc = res(t, uh, v)
  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromTFEOp{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  # Residual
  res = get_res(odeop.tfeop)
  dc = res(t, uh, v)

  # Mass
  order = get_order(odeop)
  mass = get_forms(odeop.tfeop)[1]
  ∂tNuh = ∂t(uh, Val(order))
  dc = dc + mass(t, uh, ∂tNuh, v)

  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromTFEOp{<:AbstractSemilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  # Residual
  res = get_res(odeop.tfeop)
  dc = res(t, uh, v)

  # Mass
  order = get_order(odeop)
  mass = get_forms(odeop.tfeop)[1]
  ∂tNuh = ∂t(uh, Val(order))
  dc = dc + mass(t, ∂tNuh, v)

  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOpFromTFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  # Residual
  res = get_res(odeop.tfeop)
  # Need a negative sign here:
  # residual(t, u, v) = ∑_{0 ≤ k ≤ N} form_k(t, ∂t^k[u], v) - res(t, v)
  dc = (-1) * res(t, uh, v)

  # Forms
  order = get_order(odeop)
  forms = get_forms(odeop.tfeop)
  ∂tkuh = uh
  for k in 0:order
    form = forms[k+1]
    dc = dc + form(t, ∂tkuh, v)
    if k < order
      ∂tkuh = ∂t(∂tkuh)
    end
  end

  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Algebra.allocate_jacobian(
  odeop::ODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  jacs = get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:get_order(odeop.tfeop)
    jac = jacs[k+1]
    dc = dc + jac(t, uh, du, v)
  end
  matdata = collect_cell_matrix(Ut, V, dc)
  allocate_matrix(assembler, matdata)
end

function jacobian_add!(
  J::AbstractMatrix, odeop::ODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  jacs = get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:get_order(odeop)
    w = ws[k+1]
    iszero(w) && continue
    jac = jacs[k+1]
    dc = dc + w * jac(t, uh, du, v)
  end

  if num_domains(dc) > 0
    matdata = collect_cell_matrix(Ut, V, dc)
    assemble_matrix_add!(J, assembler, matdata)
  end

  J
end

function jacobian_add!(
  J::AbstractMatrix, odeop::ODEOpFromTFEOp{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  order = get_order(odeop)
  jacs = get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:order-1
    w = ws[k+1]
    iszero(w) && continue
    jac = jacs[k+1]
    dc = dc + w * jac(t, uh, du, v)
  end

  # Special case for the mass matrix
  k = order
  w = ws[k+1]
  if !iszero(w)
    if is_form_constant(odeop, k)
      axpy_entries!(w, odeopcache.const_forms[1], J)
    else
      jac = jacs[k+1]
      dc = dc + w * jac(t, uh, du, v)
    end
  end

  if num_domains(dc) > 0
    matdata = collect_cell_matrix(Ut, V, dc)
    assemble_matrix_add!(J, assembler, matdata)
  end

  J
end

function jacobian_add!(
  J::AbstractMatrix, odeop::ODEOpFromTFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = get_assembler(odeop.tfeop)

  jacs = get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:get_order(odeop)
    w = ws[k+1]
    iszero(w) && continue
    if is_form_constant(odeop, k)
      axpy_entries!(w, odeopcache.const_forms[k+1], J)
    else
      jac = jacs[k+1]
      dc = dc + w * jac(t, uh, du, v)
    end
  end

  if num_domains(dc) > 0
    matdata = collect_cell_matrix(Ut, V, dc)
    assemble_matrix_add!(J, assembler, matdata)
  end

  J
end

#########
# Utils #
#########
# NOTE it seems that EvaluationFunction could be replaced by FEFunction. There
# is only a difference between the two functions when the underlying FESpace
# is zero mean (EvaluationFunction does not constrain the DOFs)
function _make_uh_from_us(odeop, us, Us)
  u = EvaluationFunction(Us[1], us[1])
  dus = ()
  for k in 1:get_order(odeop)
    dus = (dus..., EvaluationFunction(Us[k+1], us[k+1]))
  end
  TransientCellField(u, dus)
end
