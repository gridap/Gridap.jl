
"""
    abstract type Pushforward <: Map end

Represents a pushforward map ``F_*``, defined as ``F_*`` : V̂ -> V where
  - V̂ is a function space on the reference cell K̂ and
  - V is a function space on the physical cell K.
"""
abstract type Pushforward <: Map end

"""
    Pushforward(::ReferenceFEName, conf::Conformity)
    Pushforward(::Type{<:ReferenceFEName}, conf::Conformity)

Return the pushforward to use to map the shape functions of the given
element with conformity `conf`. For L2 conformity, the trivial
[`IdentityPiolaMap`](@ref) is always returned.

For new `ReferenceFEName`, the default pushforward is `IdentityPiolaMap`. Types
that want to change the default may only overload
[`Pushforward(::Type{<:ReferenceFEName})`](@ref).
"""
Pushforward(name::ReferenceFEName, conf::Conformity) = Pushforward(typeof(name),conf)

Pushforward(::Type{<:ReferenceFEName}, ::L2Conformity) = IdentityPiolaMap()
Pushforward(T::Type{<:ReferenceFEName}, ::Conformity) = Pushforward(T)
Pushforward(::Type{<:ReferenceFEName}) = IdentityPiolaMap()

function Arrays.lazy_map(
  k::Pushforward, ref_cell_fields::AbstractArray, pf_args::AbstractArray...
)
  lazy_map(Broadcasting(Operation(k)), ref_cell_fields, pf_args...)
end

function Arrays.evaluate!(
  cache, k::Pushforward, v_ref::Number, args...
)
  @abstractmethod
end

function evaluate!(
  cache, k::Pushforward, f_ref::AbstractVector{<:Field}, args...
)
  Broadcasting(Operation(k))(f_ref,args...)
end

function evaluate!(
  cache, k::Pushforward, f_ref::Field, args...
)
  Operation(k)(f_ref,args...)
end

function Arrays.lazy_map(
  ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{<:Pushforward}}}}
)
  cell_ref_fields, args... = a.args
  cell_ref_gradient = lazy_map(Broadcasting(∇),cell_ref_fields)
  return lazy_map(a.maps.value,cell_ref_gradient,args...)
end

function Arrays.evaluate!(
  cache,
  ::Broadcasting{typeof(gradient)},
  a::Fields.BroadcastOpFieldArray{<:Pushforward}
)
  v, pf_args... = a.args
  grad_v = Broadcasting(∇)(v)
  Broadcasting(Operation(a.op))(grad_v,pf_args...)
end

# InversePushforward

"""
    const InversePushforward{PF} = InverseMap{PF} where PF <: Pushforward

Represents the inverse of a pushforward map ``F_*``, defined as
  (``F_*``)⁻¹ : V -> V̂ where
  - V̂ is a function space on the reference cell K̂ and
  - V is a function space on the physical cell K.
"""
const InversePushforward{PF} = InverseMap{PF} where PF <: Pushforward

function Arrays.lazy_map(
  k::InversePushforward, phys_cell_fields::AbstractArray, pf_args::AbstractArray...
)
  lazy_map(Broadcasting(Operation(k)), phys_cell_fields, pf_args...)
end

function evaluate!(
  cache, k::InversePushforward, f_phys::AbstractVector{<:Field}, args...
)
  Broadcasting(Operation(k))(f_phys,args...)
end

function evaluate!(
  cache, k::InversePushforward, f_phys::Field, args...
)
  Operation(k)(f_phys,args...)
end

# Pullback

"""
    struct Pullback{PF <: Pushforward} <: Map end

Represents a pullback map ``F^*``, defined as
  ``F^*`` : V* -> V̂* where
  - V̂* is a dof space on the reference cell K̂ and
  - V* is a dof space on the physical cell K.
Its action on physical dofs σ : V -> R is defined in terms of the pushforward map ``F_*`` as:\\
σ̂ = ``F^*``(σ) := σ∘``F_*`` : V̂ -> R
"""
struct Pullback{PF <: Pushforward} <: Map
  pushforward::PF
end

function Arrays.lazy_map(
  ::typeof(evaluate),k::LazyArray{<:Fill{<:Pullback}},ref_cell_fields::AbstractArray
)
  pb = k.maps.value
  phys_cell_dofs, pf_args... = k.args
  phys_cell_fields = lazy_map(pb.pushforward,ref_cell_fields,pf_args...)
  return lazy_map(evaluate,phys_cell_dofs,phys_cell_fields)
end

function evaluate!(
  cache, k::Pullback, σ_phys::AbstractVector{<:Dof}, args...
)
  return MappedDofBasis(k.pushforward,σ_phys,args...)
end

# InversePullback

"""
    struct InversePullback{PF <: Pushforward} <: Map end

Represents the inverse of the pullback map ``(F^*)``⁻¹, defined as
  ``(F^*)``⁻¹ : V̂* -> V*
where
  - V̂* is a dof space on the reference cell K̂ and
  - V* is a dof space on the physical cell K.
Its action on reference dofs σ̂ : V̂ -> R is defined in terms of the pushforward map ``F_*`` as:\\
σ = ``(F^*)``⁻¹(σ̂) := σ̂∘``(F_*)``⁻¹ : V -> R
"""
const InversePullback{PB} = InverseMap{PB} where PB <: Pullback

function Arrays.lazy_map(
  ::typeof(evaluate), k::LazyArray{<:Fill{<:InversePullback}}, phys_cell_fields::AbstractArray
)
  pb = inverse_map(k.maps.value)
  ref_cell_dofs, pf_args... = k.args
  ref_cell_fields = lazy_map(inverse_map(pb.pushforward), phys_cell_fields, pf_args...)
  return lazy_map(evaluate,ref_cell_dofs,ref_cell_fields)
end

function evaluate!(
  cache, k::InversePullback, σ_ref::AbstractVector{<:Dof}, args...
)
  pb = inverse_map(k)
  return MappedDofBasis(inverse_map(pb.pushforward),σ_ref,args...)
end

##############
# Piola maps #
##############

# In what follows,
# - F is the geometrical map F:K̂->K
# - Jt = ∇F = (Jac(F))ᵀ

"""
    struct IdentityPiolaMap <: Pushforward
"""
struct IdentityPiolaMap <: Pushforward end # φ̂ -> φ = φ̂∘F⁻¹

# ContraVariantPiolaMap

"""
    struct ContraVariantPiolaMap <: Pushforward
"""
struct ContraVariantPiolaMap <: Pushforward end

function evaluate!( # φ̂ -> φ = (|det(J)|⁻¹J φ̂)∘F⁻¹
  cache, ::ContraVariantPiolaMap, v_ref::Number, Jt::Number
)
  idetJ = 1. / meas(Jt)
  return v_ref ⋅ (idetJ * Jt)
end

function evaluate!( # φ -> φ̂ = |det(J)| J⁻¹ φ∘F
  cache, ::InversePushforward{ContraVariantPiolaMap}, v_phys::Number, Jt::Number
)
  detJ = meas(Jt)
  return v_phys ⋅ (detJ * pinvJt(Jt))
end

# TODO: Should this be here? Probably not...

function Fields.DIV(f::LazyArray{<:Fill})
  df = Fields.DIV(f.args[1])
  k  = f.maps.value
  lazy_map(k,df)
end

function Fields.DIV(a::LazyArray{<:Fill{typeof(linear_combination)}})
  i_to_basis  = Fields.DIV(a.args[2])
  i_to_values = a.args[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

function Fields.DIV(f::LazyArray{<:Fill{Broadcasting{Operation{ContraVariantPiolaMap}}}})
  ϕrgₖ = f.args[1]
  return lazy_map(Broadcasting(divergence),ϕrgₖ)
end

function Fields.DIV(f::Fill{<:Fields.BroadcastOpFieldArray{ContraVariantPiolaMap}})
  ϕrgₖ = f.value.args[1]
  return Fill(Broadcasting(divergence)(ϕrgₖ),length(f))
end

# CoVariantPiolaMap

"""
    struct CoVariantPiolaMap <: Pushforward
"""
struct CoVariantPiolaMap <: Pushforward end

function evaluate!( # φ̂ -> φ = (J⁻ᵀ φ̂)∘F⁻¹
  cache, ::CoVariantPiolaMap, v_ref::Number, Jt::Number
)
  return v_ref ⋅ transpose(pinvJt(Jt))
end

function evaluate!( # φ -> φ̂ = Jᵀ φ∘F
  cache, ::InversePushforward{CoVariantPiolaMap}, v_phys::Number, Jt::Number
)
  return v_phys ⋅ transpose(Jt)
end

# CoContraVariantPiolaMap

"""
    struct CoContraVariantPiolaMap <: Pushforward

The mixed Piola map `φ̂ ↦ φ = det(J)⁻¹ J⁻ᵀ φ̂ Jᵀ`, covariant on the first index
and contravariant on the second, for matrix-valued fields whose
*normal-tangential* components are the continuous ones. It is the map of the
Gopalakrishnan--Lederer--Schöberl element.
Does not preserve symmetry, but preserves the trace.
"""
struct CoContraVariantPiolaMap <: Pushforward end

function evaluate!( # φ̂ -> φ = (det(J)⁻¹ J⁻ᵀ φ̂ Jᵀ)∘F⁻¹
  cache, ::CoContraVariantPiolaMap, v_ref::Number, Jt::Number
)
  return ((1. / det(Jt)) * pinvJt(Jt)) ⋅ v_ref ⋅ Jt
end

function evaluate!( # φ -> φ̂ = det(J) Jᵀ φ∘F J⁻ᵀ
  cache, ::InversePushforward{CoContraVariantPiolaMap}, v_phys::Number, Jt::Number
)
  return (det(Jt) * Jt) ⋅ v_phys ⋅ pinvJt(Jt)
end

# DoubleContraVariantPiolaMap

struct DoubleContraVariantPiolaMap <: Pushforward end

function evaluate!( # φ̂ -> φ = (det(J)⁻² J φ̂ Jᵀ)∘F⁻¹
  cache, ::DoubleContraVariantPiolaMap, v_ref::Number, Jt::Number
)
  _Jt = (1. / det(Jt)) * Jt
  return congruent_prod(v_ref, _Jt) # symmetry stable _Jtᵀ ⋅ v_ref ⋅ _Jt
end

function evaluate!( # φ -> φ̂ = det(J)² J⁻¹ φ∘F J⁻ᵀ
  cache, ::InversePushforward{DoubleContraVariantPiolaMap}, v_phys::Number, Jt::Number
)
  iJt = det(Jt) * pinvJt(Jt)
  return congruent_prod(v_phys, iJt) # symmetry stable iJtᵀ ⋅ v_ref ⋅ iJt
end

function Fields.DIV(f::LazyArray{<:Fill{Broadcasting{Operation{DoubleContraVariantPiolaMap}}}})
  ϕrgₖ, Jt = f.args
  return lazy_map(ContraVariantPiolaMap(), lazy_map(Broadcasting(divergence), ϕrgₖ), Jt)
end

function Fields.DIV(f::Fill{<:Fields.BroadcastOpFieldArray{DoubleContraVariantPiolaMap}})
  ϕrgₖ, Jt = f.value.args
  divϕ = Broadcasting(Operation(ContraVariantPiolaMap()))(Broadcasting(divergence)(ϕrgₖ), Jt)
  return Fill(divϕ, length(f))
end

# DoubleCoVariantPiolaMap

struct DoubleCoVariantPiolaMap <: Pushforward end

function evaluate!( # φ̂ -> φ = (J⁻ᵀ φ̂ J⁻¹)∘F⁻¹
  cache, ::DoubleCoVariantPiolaMap, v_ref::Number, Jt::Number
)
  iJt = pinvJt(Jt)
  return congruent_prod(v_ref, transpose(iJt)) # symmetry stable iJt ⋅ v_ref ⋅ iJtᵀ
end

function evaluate!( # φ -> φ̂ = Jᵀ φ∘F J
  cache, ::InversePushforward{DoubleCoVariantPiolaMap}, v_phys::Number, Jt::Number
)
  return congruent_prod(v_phys, transpose(Jt)) # symmetry stable Jt ⋅ v_phys ⋅ Jtᵀ
end


###################
# DOF scaling API #
###################

"""
    get_dofscale_setter_function(reffe::ReferenceFE, pushforward::Pushforward)

Return a function `(dofscale, face_own_dofs, face_meshsize) -> _` that
sets in place in `dofscale` the scaling with `h` that the DOF `dof` of `reffe`
undergoes when mapped to the physical space, where `h` is a local meshsize
estimate of the face owning `dof`.

Returned function arguments:
- `dofscale`: vector of floats of length `num_dofs(reffe)`,
- `face_own_dofs`: result of [`get_face_own_dofs(reffe)`](@ref get_face_own_dofs),
- `face_meshsize`: vector of floats of length `num_faces(reffe)`, its `i`th entry is the local meshize estimate for the DOFs in `face_own_dofs[i]` (if non-empty).

By default, the scaling function is deduced from the `reffe`'s `Pushforward` (Piola map).

For example, the standard div-conforming FEs (Raviart-Thomas, BDM) map all their
DOFs with the contravariant Piola map which applies `J -> |det(J)|⁻¹ J`, that
scales like `h^D * (1/h) = h^(D-1)`. So `get_dofscale_setter_function` returns
something equivalent to

    @inline function scale_setter!(dofscale, face_own_dofs, face_meshsize)
      for (face_dofs, h) in zip(face_own_dofs, face_meshsize)
        for dof in face_dofs
          dofscale[dof] = h^(D-1)
        end
      end
    end

This function is responsible for the correct DOF scaling if using the `scale_dof`
kwarg of the [`FESpace`](@ref) constructor.
"""
function get_dofscale_setter_function(::ReferenceFE{D}, pushforward::Pushforward) where D
  scaler = _scaling_function(pushforward, D)

  # avoid boxing scaler in the closure by making it constant
  # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
  let scaler=scaler

    # The returned "scale_setter!" anonymous funtion:
    @inline function(dofscale, face_own_dofs, face_meshsize)
      for (face_dofs, h) in zip(face_own_dofs, face_meshsize)
        hp = scaler(h)
        for dof in face_dofs
          dofscale[dof] = hp
        end
      end
    end

  end
end

_scaling_function(::Pushforward, D::Int) = @abstractmethod
_scaling_function(::IdentityPiolaMap, D::Int)  = h -> 1
_scaling_function(::CoVariantPiolaMap, D::Int) = h -> h
_scaling_function(::DoubleCoVariantPiolaMap, D::Int) = h -> h^2
_scaling_function(::ContraVariantPiolaMap, D::Int) = let D=D
  h -> h^(D-1)
end
_scaling_function(::DoubleContraVariantPiolaMap, D::Int) = let D=D
  h -> h^(2D-2)
end
_scaling_function(::CoContraVariantPiolaMap, D::Int) = let D=D
  h -> h^D
end

# The scale setter of a reffe whose DoF `i` scales like `h^exponent[i]`
function _dofscale_setter_from_exponents(exponent::AbstractVector{<:Integer})
  let exponent=exponent
    @inline function(dofscale, face_own_dofs, face_meshsize)
      for (face_dofs, h) in zip(face_own_dofs, face_meshsize)
        e, he = 0, one(h)
        for dof in face_dofs
          if exponent[dof] != e
            e = exponent[dof]
            he = h^e
          end
          dofscale[dof] = he
        end
      end
    end
  end
end

# The setter of an element whose every DoF is a moment over the face owning it,
# scaling like the `d`-measure of its `d`-face: `h^d`
function _face_dim_dofscale_setter(reffe::ReferenceFE)
  exponent = zeros(Int, num_dofs(reffe))
  for (face_dofs, d) in zip(get_face_own_dofs(reffe), get_facedims(get_polytope(reffe)))
    exponent[face_dofs] .= d
  end
  _dofscale_setter_from_exponents(exponent)
end

################################################################################
# Change of basis of the edge-scaled elements (HHJ, Regge, GLS)
#
# The cell-local map. The mesh-level `compute_cell_bases_changes` methods live in
# src/FESpaces/Pullbacks.jl.

# The change of basis of an element whose edge DoFs are *preserved up to one
# positive scalar per edge*: diagonal, with `‖J t̂ₑ‖^{±1}` on the DoFs of edge `e`,
# times `(-1)ⁱ` when the cell traverses that edge against the global direction, and
# `1` on the interior DoFs. `transposed_inverse` selects `P⁻ᵀ` over `P`.
#
# This covers every element whose edge DoF has the form `∫ₑ (a⋅M⋅b) μᵢ ds` for a
# pair of directions carried dually by that element's push-forward — a normal and a
# normal under the double contravariant map, two tangents under the double
# covariant one, one of each under the co-contravariant one. In every case the
# Jacobians cancel completely,
#
#     a⋅M⋅b = (â⋅M̂⋅b̂) / ‖J t̂ₑ‖²,   ds = ‖J t̂ₑ‖ dŝ   ⟹   F∗(ℓ^{e,i}) = ℓ̂^{e,i} / ‖J t̂ₑ‖.
#
# That is not a coincidence. A Piola map is *chosen* so that its element's DoF is
# invariant, so what is left over cannot depend on which pairing was picked — only
# on the fact that the DoF is an edge moment against a degree-`i` weight. The
# elements still need their own `compute_cell_bases_changes`, which dispatches on
# the reference FE name and the push-forward.
#
# Since `a⋅M⋅b` is quadratic in the directions, or bilinear with both of them
# flipping, none of these elements needs a normal or tangent sign convention; the
# only orientation effect is the parity of the weight, which is the `(-1)ⁱ` above.
# At order 0 even that disappears.
struct EdgeScalingChangeOfBasis <: Map
  tangents::Vector{VectorValue{2,Float64}}
  edge_dofs::Vector{Vector{Int}}
  ndofs::Int
  transposed_inverse::Bool
end

function EdgeScalingChangeOfBasis(reffe::ReferenceFE, transposed_inverse::Bool)
  p = get_polytope(reffe)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  edge_dofs = [own[nv+e] for e in 1:num_faces(p, 1)]
  return EdgeScalingChangeOfBasis(
    get_edge_tangent(p), edge_dofs, num_dofs(reffe), transposed_inverse
  )
end

function return_cache(k::EdgeScalingChangeOfBasis, Jt, σ)
  CachedArray(zeros(Float64, k.ndofs, k.ndofs))
end

function evaluate!(cache, k::EdgeScalingChangeOfBasis, Jt, σ)
  setsize!(cache, (k.ndofs, k.ndofs))
  M = cache.array
  fill!(M, zero(eltype(M)))
  for i in 1:k.ndofs
    M[i, i] = 1.0
  end

  for e in eachindex(k.edge_dofs)
    L = norm(k.tangents[e] ⋅ Jt)   # ‖J t̂ₑ‖, the edge length ratio
    reversed = σ[e] < 0
    for (i, d) in enumerate(k.edge_dofs[e])
      # the i-th DoF of the edge carries the weight of degree i-1
      s = ifelse(reversed && isodd(i - 1), -1.0, 1.0)
      M[d, d] = ifelse(k.transposed_inverse, s / L, s * L)
    end
  end

  return M
end
