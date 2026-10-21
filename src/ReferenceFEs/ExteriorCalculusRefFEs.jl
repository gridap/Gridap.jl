############################################
# (Proxied) Form valued FEEC reference FEs #
############################################

function ReferenceFE(p::Polytope{D},F::Symbol,r,k,T::Type;
  rotate_90=false, DG_calc=true, nodal::Bool=false, kwargs...) where D

  _DG_calc = DG_calc || T <: ExteriorFormValue
  cart_prod = !_DG_calc && k in (0,D) && _is_cart_prod(T)
  T <: ExteriorFormValue && @check _form_degree_dim(T) == (k,D) """
    The value type $T is a $(first(_form_degree_dim(T)))-form in dimension \
    $(last(_form_degree_dim(T))), but a $k-form in dimension $D was asked for.
    """
  V = T <: ExteriorFormValue ? T : (_DG_calc ? ExteriorFormValue{k,D,T} : T)
  FEEC_space_definition_checks(Val(D),V,r,k,F,rotate_90; cart_prod, element=true)

  # what FEEC_space_definition_checks cannot check, one involves the polytope,
  # the other the choice of DOF functionals.
  @check ((F in (:P,:P⁻)) ? is_simplex(p) : is_n_cube(p)) """
    The $F family lives on $(F in (:P,:P⁻) ? "simplices" : "n-cubes"), got $p
    """
  nodal && @check k in (0,D) "Nodal DOFs are only meaningful for 0-forms and D-forms, got k=$k, D=$D"

  VV = _DG_calc ? make_concretetype(V) : V

  if F == :S && k == D # 𝓢ᵣΛᴰ ≡ 𝓟ᵣ space, which no modal name carries
    @notimplementedif !nodal "Modal SrΛᴰ elements are not implemented, please use nodal ones (pass `nodal=true`)."
    return ReferenceFE(p, lagrangian, VV, r; space=:P, kwargs...)
  end

  (name, order) = FEEC_name_and_order(F, r, k, Val(D), nodal, rotate_90)
  # `name` is a legacy (but cell agnostic) `ReferenceFEName`, like lagrangian,
  # nedelec, etc. So this dispatches to the file of `name`, but for
  # V<:ExteriorFormValue, the FEEC_reference_fe defined below will be called next.
  ReferenceFE(p, name, VV, order; kwargs...)
end


###########################################################
# Form valued (ExteriorFormValue) FEEC reference elements #
###########################################################

"""
    FEEC_reference_fe(name, ::Type{V<:ExteriorFormValue}, p::Polytope{D}, order; kwargs...)

Build the form valued FEEC reference element that `name` denotes, with shape
functions valued in the [`ExteriorFormValue`](@ref) type `V`.

This is what the `*RefFE` constructors below forward to when handed a form
value type; it is not called directly.
"""
function FEEC_reference_fe(
  name::ReferenceFEName, ::Type{V}, p::Polytope{D}, order::Integer;
  change_dof=true, poly_type=_mom_reffe_default_PT(p), mom_poly_type=poly_type
) where {V<:ExteriorFormValue,D}

  VV = make_concretetype(V)
  k, Dv = _form_degree_dim(VV)
  T = eltype(VV)

  @check Dv == D "The value type $V lives in dimension $Dv, but the polytope is $D dimensional"
  @check 0 ≤ k ≤ D "The form degree must verify 0 ≤ k ≤ D, got k=$k, D=$D"
  @check k in FEEC_form_degrees(name, D) """
    $name names a space of $(join(FEEC_form_degrees(name, D), " or "))-forms \
    in dimension $D, but the value type $V is a $k-form. \
    In dimension $D the $k-forms are named by $(_FEEC_degree_names(k, D)).
    """
  F, r = FEEC_family_and_order(name, p, order, k)

  @check (F in (:P,:P⁻) && is_simplex(p)) || (F == :Q⁻ && is_n_cube(p)) """
    Form valued *moment based* reference FEs are available for the P and P⁻
    families on simplices, and for the Q⁻ family on n-cubes; got the $F family
    on $p.
    """

  PT, MPT = poly_type, mom_poly_type
  # the value type asks for the form valued bases, no flag needed
  prebasis = FEEC_poly_basis(Val(D), VV, r, k, F, PT)

  # A d-face carries dofs for d ≥ k
  moments = Tuple[]
  for d in k:D
    mom_r, mom_F = FEEC_face_mom_order_and_family(Val(F), r, k, d)
    isnothing(mom_F) && continue
    mom_k = d-k
    μV = ExteriorFormValue{mom_k,d,T}
    μ = FEEC_poly_basis(Val(d), μV, mom_r, mom_k, mom_F, MPT)
    iszero(length(μ)) && continue
    push!(moments, (get_dimrange(p,d), form_moment, μ))
  end

  conf = PΛConformity()
  change_dof = _validate_change_dof(change_dof, prebasis, p, conf)
  MomentBasedReferenceFE(name, p, prebasis, moments, conf; change_dof)
end

get_form_degree(reffe::ReferenceFE) = get_form_degree(return_type(get_shapefuns(reffe)))

"""
    form_moment(φ, μ, ds) -> Field-like

The AFW moment integrand ``φ → ∫( ι^*φ ∧ μ )df``, where `ι` is the inclusion
from the cell to the face `f` of `ds`.
"""
function form_moment(φ, μ, ds::FaceMeasure)
  # the trace of φ on the face replaces the extention of μ to the cell.
  Jt  = Broadcasting(∂)(ds.fmaps[ds.face])
  ιφ  = Broadcasting(Operation(pullback))(φ, Jt)
  wdg = Broadcasting(Operation(∧))(ιφ, μ)
  Broadcasting(Operation(flat_vol_coeff))(wdg)
end

# ─────────────────────────────────────────────────────────────────────────────
# helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    FEEC_name_and_order(F::Symbol, r, k, ::Val{D}, nodal, rotate_90) -> (name, order)

The `(F,r,k) ↦ (name, order)` table of the Periodic Table of the Finite
Elements, shared by the form valued and the vector proxied elements -- the two
agree on every cell, `rotate_90` only picking between the two names a 2D space
has.

This logic should keep consistent with the Table in the ReferenceFEs
documentation. `FEEC_family_and_order` below is its inverse, name by
name, up to `𝓟ᵣΛ⁰ = 𝓟ᵣ⁻Λ⁰` and `𝓢ᵣΛᴰ = 𝓟ᵣΛᴰ`, where two families name one space.
"""
function FEEC_name_and_order(F::Symbol, r, k, ::Val{D}, nodal::Bool, rotate_90::Bool) where D
  # The single component degrees are the scalar ends of the table, nodal or
  # moment based, and are tested first so that D=1 lands on them.
  if k in (0,D)
    scalar = F == :S ? serendipity : lagrangian
    # `𝓟ᵣ⁽⁻⁾Λ⁰ = 𝓟ᵣ`, `𝓟ᵣ⁻Λᴰ = 𝓟ᵣ₋₁`, `𝓟ᵣΛᴰ = 𝓟ᵣ`, and likewise for 𝓠/𝓢
    order = (iszero(k) || F in (:P,:S)) ? r : r-1
    return (nodal ? scalar : ModalScalar(scalar), order)
  elseif k == 1 && !rotate_90
    F == :S && @notimplemented "Serendipity Nedelec is not implemented yet"
    return F == :P ? (nedelec2, r) : (nedelec1, r-1)
  elseif k == D-1 # in 2D, the k=1 space reached with the `rotate_90` proxy
    F == :S && @notimplemented "BDM on n-cubes is not implemented yet"
    return F == :P ? (bdm, r) : (raviart_thomas, r-1)
  end
  @unreachable "No reference FE name stands for the $k-forms in dimension $D"
end

# `F_sΛ^j` is defined for s ≥ 0, except that the trimmed family starts at s = 1.
# The threshold is the one of the *moment* family returned, not of `F`: a `:P`
# element takes trimmed `:P⁻` moments, and a `:P⁻` element full `:P` ones.
function FEEC_face_mom_order_and_family(::Val{:P}, r, k, d)
  mom_r = r+k-d
  mom_r ≥ 1 ? (mom_r, :P⁻) : (-1, nothing)
end
function FEEC_face_mom_order_and_family(::Val{:P⁻}, r, k, d)
  mom_r = r+k-d-1
  mom_r ≥ 0 ? (mom_r, :P)  : (-1, nothing)
end
function FEEC_face_mom_order_and_family(::Val{:Q⁻}, r, k, d)
  mom_r = r-1
  mom_r ≥ 0 ? (mom_r, :Q⁻) : (-1, nothing)
end

FEEC_form_degrees(::Nedelec, D) = (1,)
FEEC_form_degrees(::RaviartThomas, D) = (D-1,)
FEEC_form_degrees(::BDM, D) = (D-1,)
FEEC_form_degrees(::ModalScalar, D) = (0,D)

function FEEC_family_and_order(::Nedelec{1}, p, order, k)
  trimmed = is_simplex(p) ? :P⁻ : :Q⁻
  (trimmed, order+1)
end
function FEEC_family_and_order(::Nedelec{2}, p, order, k)
  @check is_simplex(p) "Nedelec reference elements of the second kind are only defined on simplices"
  (:P, order)
end
function FEEC_family_and_order(::RaviartThomas, p, order, k)
  trimmed = is_simplex(p) ? :P⁻ : :Q⁻
  (trimmed, order+1)
end
function FEEC_family_and_order(::BDM, p, order, k)
  @check is_simplex(p) "BDM reference elements are only defined on simplices"
  (:P, order)
end
function FEEC_family_and_order(::ModalScalar{Lagrangian}, p, order, k)
  trimmed = is_simplex(p) ? :P⁻ : :Q⁻
  iszero(k) ? (is_simplex(p) ? :P : trimmed, order) : (trimmed, order+1)
end
function FEEC_family_and_order(::ModalScalar{Serendipity}, p, order, k)
  @check is_n_cube(p) "Serendipity reference elements are only defined on n-cubes"
  (:S, order)
end

# For the error message of the form degree check above.
function _FEEC_degree_names(k, D)
  k in (0,D)  && return "`modal_lagrangian`, `lagrangian` or `serendipity`"
  k == 1      && return D == 2 ? "`nedelec`, `nedelec2`, `raviart_thomas` or `bdm`" :
                                 "`nedelec` or `nedelec2`"
  k == D-1    && return "`raviart_thomas` or `bdm`"
  "no reference FE name"
end

