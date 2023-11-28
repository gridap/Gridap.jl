######################
# ODEOpFromFEOpCache #
######################
"""
    struct ODEOpFromFEOpCache <: GridapType

Structure that stores the `TransientFESpace` and cache of a
`TransientFEOperator`, as well as the jacobian matrices and forcing term if
they are constant.
"""
mutable struct ODEOpFromFEOpCache <: GridapType
  Us
  Uts
  fe_cache
  jacs
  jacvec
  forcing
end

#################
# ODEOpFromFEOp #
#################
"""
    struct ODEOpFromFEOp <: ODEOperator end

Wrapper that transforms a `TransientFEOperator` into an `ODEOperator`, i.e.
takes `residual(t, uh, ∂t(uh), ..., ∂t^N(uh), vh)` and returns
`residual(t, uf, ∂t(uf), ..., ∂t^N(uf))`, where `uf` represent the free values
of the `EvaluationFunction` `uh`.
"""
struct ODEOpFromFEOp{C} <: ODEOperator{C}
  fe_op::TransientFEOperator{C}
end

# ODEOperator interface
Polynomials.get_order(op::ODEOpFromFEOp) = get_order(op.fe_op)

function allocate_cache(
  op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}}
)
  order = get_order(op)
  is_masslinear = ODEOperatorType(op) <: AbstractMassLinearODE

  Ut = get_trial(op.fe_op)
  U = allocate_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for k in 1:order
    Uts = (Uts..., ∂t(Uts[k]))
    Us = (Us..., allocate_space(Uts[k+1]))
  end
  fe_cache = allocate_cache(op.fe_op)

  V = get_test(op.fe_op)
  v = get_fe_basis(V)
  uh = _make_uh_from_us(op, us, Us)

  jacs = ()
  use_jacvec = false
  for k in 0:order
    jac = nothing
    if is_jacobian_constant(op, k)
      matdata = _matdata_jacobian(op.fe_op, t, uh, k, 1)
      jac = assemble_matrix(get_assembler(op.fe_op), matdata)

      use_jacvec = true
    end
    jacs = (jacs..., jac)
  end

  jacvec = nothing
  if is_masslinear && use_jacvec
    res = get_res(op.fe_op)
    vecdata = collect_cell_vector(V, res(t, uh, v))
    jacvec = allocate_vector(get_assembler(op.fe_op), vecdata)
  end

  forcing = nothing
  if is_masslinear && is_forcing_constant(op)
    res = get_res(op.fe_op)
    vecdata = collect_cell_vector(V, res(t, uh, v))
    forcing = assemble_vector(get_assembler(op.fe_op), vecdata)
  end

  ODEOpFromFEOpCache(Us, Uts, fe_cache, jacs, jacvec, forcing)
end

function update_cache!(cache, op::ODEOpFromFEOp, t::Real)
  Us = ()
  for k in 0:get_order(op)
    Us = (Us..., evaluate!(cache.Us[k+1], cache.Uts[k+1], t))
  end
  cache.Us = Us
  cache.fe_cache = update_cache!(cache.fe_cache, op.fe_op, t)
  cache
end

############
# Residual #
############
function Algebra.allocate_residual(
  op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  uh = _make_uh_from_us(op, us, cache.Us)

  V = get_test(op.fe_op)
  v = get_fe_basis(V)
  res = get_res(op.fe_op)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  allocate_vector(get_assembler(op.fe_op), vecdata)
end

function Algebra.residual!(
  r::AbstractVector, op::ODEOpFromFEOp{NonlinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  uh = _make_uh_from_us(op, us, cache.Us)

  V = get_test(op.fe_op)
  v = get_fe_basis(V)
  res = get_res(op.fe_op)
  vecdata = collect_cell_vector(V, res(t, uh, v))
  assemble_vector!(r, get_assembler(op.fe_op), vecdata)
  r
end

function Algebra.residual!(
  r::AbstractVector, op::ODEOpFromFEOp{MassLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache; include_highest::Bool=true
)
  mass_stored = !isnothing(cache.jacs[end])
  forcing_stored = !isnothing(cache.forcing)
  if !(mass_stored && forcing_stored)
    V = get_test(op.fe_op)
    v = get_fe_basis(V)
    uh = _make_uh_from_us(op, us, cache.Us)
  end

  fill!(r, zero(eltype(r)))

  if include_highest
    if mass_stored
      mul!(cache.jacvec, cache.jacs[end], us[end])
      axpy!(1, cache.jacvec, r)
    else
      mass = get_mass(op.fe_op)
      vecdata = collect_cell_vector(V, mass(t, uh, v))
      assemble_vector_add!(r, get_assembler(op.fe_op), vecdata)
    end
  end

  if forcing_stored
    axpy!(1, cache.forcing, r)
  else
    res = get_res(op.fe_op)
    vecdata = collect_cell_vector(V, res(t, uh, v))
    assemble_vector_add!(r, get_assembler(op.fe_op), vecdata)
  end

  r
end

function Algebra.residual!(
  r::AbstractVector, op::ODEOpFromFEOp{LinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache; include_highest::Bool=true
)
  forms_stored = !any(isnothing, cache.jacs)
  forcing_stored = !isnothing(cache.forcing)
  if !(forms_stored && forcing_stored)
    V = get_test(op.fe_op)
    v = get_fe_basis(V)
    uh = _make_uh_from_us(op, us, cache.Us)
  end

  fill!(r, zero(eltype(r)))

  order = get_order(op)
  order_max = include_highest ? order : order - 1
  for k in 0:order_max
    if !isnothing(cache.jacs[k+1])
      mul!(cache.jacvec, cache.jacs[k+1], us[k+1])
      axpy!(1, cache.jacvec, r)
    else
      form = get_forms(op.fe_op)[order+1-k]
      vecdata = collect_cell_vector(V, form(t, uh, v))
      assemble_vector_add!(r, get_assembler(op.fe_op), vecdata)
    end
  end

  if forcing_stored
    axpy!(1, cache.forcing, r)
  else
    res = get_res(op.fe_op)
    vecdata = collect_cell_vector(V, res(t, uh, v))
    assemble_vector_add!(r, get_assembler(op.fe_op), vecdata)
  end

  r
end

############
# Jacobian #
############
function Algebra.allocate_jacobian(
  op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  uh = _make_uh_from_us(op, us, cache.Us)

  _matdata = ()
  for k in 0:get_order(op.fe_op)
    _matdata = (_matdata..., _matdata_jacobian(op.fe_op, t, uh, k, 0))
  end

  matdata = _vcat_matdata(_matdata)
  allocate_matrix(get_assembler(op.fe_op), matdata)
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  cache
)
  if !isnothing(cache.jacs[k+1])
    axpy_sparse!(γ, cache.jacs[k+1], J)
  else
    uh = _make_uh_from_us(op, us, cache.Us)
    matdata = _matdata_jacobian(op.fe_op, t, uh, k, γ)
    assemble_matrix_add!(J, get_assembler(op.fe_op), matdata)
  end

  J
end

function jacobians!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  need_assembly = any(isnothing, cache.jacs)
  if need_assembly
    uh = _make_uh_from_us(op, us, cache.Us)
  end

  _matdata = ()
  for k in 0:get_order(op)
    γ = γs[k+1]
    if !iszero(γ)
      if !isnothing(cache.jacs[k+1])
        axpy_sparse!(γ, cache.jacs[k+1], J)
      else
        _matdata = (_matdata..., _matdata_jacobian(op.fe_op, t, uh, k, γ))
      end
    end
  end

  if need_assembly
    matdata = _vcat_matdata(_matdata)
    assemble_matrix_add!(J, get_assembler(op.fe_op), matdata)
  end

  J
end

############
# Constant #
############
function is_jacobian_constant(op::ODEOpFromFEOp, k::Integer)
  is_jacobian_constant(op.fe_op, k)
end

function is_forcing_constant(op::ODEOpFromFEOp{<:AbstractMassLinearODE})
  is_forcing_constant(op.fe_op)
end

#########
# Utils #
#########
function _matdata_jacobian(
  op::TransientFEOperatorTypes,
  t::Real, uh::TransientCellField,
  k::Integer, γ::Real
)
  Ut = evaluate(get_trial(op), nothing)
  V = get_test(op)
  du = get_trial_fe_basis(Ut)
  v = get_fe_basis(V)
  jac = get_jacs(op)[k+1]
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

function _make_uh_from_us(op, us, Us)
  u = EvaluationFunction(Us[1], us[1])

  dus = ()
  for k in 1:get_order(op)
    dus = (dus..., EvaluationFunction(Us[k+1], us[k+1]))
  end

  if first(Us) isa MultiFieldFESpace
    TransientCellFieldType = TransientMultiFieldCellField
  else
    TransientCellFieldType = TransientCellField
  end
  TransientCellFieldType(u, dus)
end

function axpy_sparse!(α, A, B)
  # TODO optimise the sum of sparse matrices
  # This is surprisingly better than axpy!(α, A, B)
  @. B += α * A
end
