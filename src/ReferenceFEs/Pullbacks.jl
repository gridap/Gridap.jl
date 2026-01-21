
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
  return congruent_prod(v_ref, transpose(Jt)) # symmetry stable Jt ⋅ v_phys ⋅ Jtᵀ
end
