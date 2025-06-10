
"""
    abstract type Pushforward <: Map end

Represents a pushforward map F*, defined as
  F* : V̂ -> V
where
  - V̂ is a function space on the reference cell K̂ and
  - V is a function space on the physical cell K.
"""
abstract type Pushforward <: Map end

abstract type PushforwardRefFE <: ReferenceFEName end

Pushforward(::Type{<:PushforwardRefFE}) = @abstractmethod
Pushforward(name::PushforwardRefFE) = Pushforward(typeof(name))

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

Represents the inverse of a pushforward map F*, defined as
  (F*)⁻¹ : V -> V̂
where
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

Represents a pullback map F**, defined as
  F** : V* -> V̂*
where
  - V̂* is a dof space on the reference cell K̂ and
  - V* is a dof space on the physical cell K.
Its action on physical dofs σ : V -> R is defined in terms of the pushforward map F* as
 σ̂ = F**(σ) := σ∘F* : V̂ -> R
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

Represents the inverse of the pullback map F**, defined as
  (F**)⁻¹ : V̂* -> V*
where
  - V̂* is a dof space on the reference cell K̂ and
  - V* is a dof space on the physical cell K.
Its action on reference dofs σ̂ : V̂ -> R is defined in terms of the pushforward map F* as
σ = (F**)⁻¹(σ̂) := σ̂∘(F*)⁻¹ : V -> R
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

# ContraVariantPiolaMap

"""
    struct ContraVariantPiolaMap <: Pushforward
"""
struct ContraVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::ContraVariantPiolaMap, v_ref::Number, Jt::Number
)
  idetJ = 1. / meas(Jt)
  return v_ref ⋅ (idetJ * Jt)
end

function evaluate!(
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

struct CoVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::CoVariantPiolaMap, v_ref::Number, Jt::Number
)
  return v_ref ⋅ transpose(pinvJt(Jt))
end

function evaluate!(
  cache, ::InversePushforward{CoVariantPiolaMap}, v_phys::Number, Jt::Number
)
  return v_phys ⋅ transpose(Jt)
end

# DoubleContraVariantPiolaMap

struct DoubleContraVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::DoubleContraVariantPiolaMap, v_ref::Number, Jt::Number
)
  _Jt = (1. / det(Jt)) * Jt
  res = transpose(_Jt) ⋅ v_ref ⋅ _Jt
  return symmetric_part(res)
end

function evaluate!(
  cache, ::InversePushforward{DoubleContraVariantPiolaMap}, v_phys::Number, Jt::Number
)
  iJt = det(Jt) * pinvJt(Jt)
  res = transpose(iJt) ⋅ v_phys ⋅ iJt
  return symmetric_part(res)
end

# DoubleCoVariantPiolaMap

struct DoubleCoVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::DoubleCoVariantPiolaMap, v_ref::Number, Jt::Number
)
  iJt = pinvJt(Jt)
  return iJt ⋅ v_ref ⋅ transpose(iJt)
end

function evaluate!(
  cache, ::InversePushforward{DoubleCoVariantPiolaMap}, v_phys::Number, Jt::Number
)
  return Jt ⋅ v_phys ⋅ transpose(Jt)
end
