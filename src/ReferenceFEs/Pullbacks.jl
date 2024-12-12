
"""
    abstract type Pushforward <: Map end

Represents a pushforward map F*, defined as 
  F* : V̂ -> V
where 
  - V̂ is a function space on the reference cell K̂ and 
  - V is a function space on the physical cell K.
"""
abstract type Pushforward <: Map end

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

function Arrays.lazy_map(
  ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{<:Pushforward}}}}
)
  cell_ref_fields, args = a.args
  cell_ref_gradient = lazy_map(Broadcasting(∇),cell_ref_fields)
  return lazy_map(a.maps.value,cell_ref_gradient,args...)
end

function Arrays.evaluate!(
  cache,
  ::Broadcasting{typeof(∇)},
  a::Fields.BroadcastOpFieldArray{<:Pushforward}
)
  v, Jt, sign_flip = a.args
  ∇v = Broadcasting(∇)(v)
  k = ContraVariantPiolaMap()
  Broadcasting(Operation(k))(∇v,Jt,sign_flip)
end

# InversePushforward

"""
    struct InversePushforward{PF <: Pushforward} <: Map end

Represents the inverse of a pushforward map F*, defined as 
  (F*)^-1 : V -> V̂
where 
  - V̂ is a function space on the reference cell K̂ and 
  - V is a function space on the physical cell K.
"""
struct InversePushforward{PF} <: Map
  pushforward::PF
  function InversePushforward(pushforward::Pushforward)
    PF = typeof(pushforward)
    new{PF}(pushforward)
  end
end

Arrays.inverse_map(pf::Pushforward) = InversePushforward(pf)
Arrays.inverse_map(ipf::InversePushforward) = ipf.pushforward

function Arrays.lazy_map(
  k::InversePushforward, phys_cell_fields::AbstractArray, pf_args::AbstractArray...
)
  lazy_map(Broadcasting(Operation(k)), phys_cell_fields, pf_args...)
end

function Arrays.return_cache(
  k::InversePushforward, v_phys::Number, args...
)
  mock_basis(::VectorValue{D,T}) where {D,T} = one(TensorValue{D,D,T})
  v_ref_basis = mock_basis(v_phys)
  pf_cache = return_cache(k.pushforward,v_ref_basis,args...)
  return v_ref_basis, pf_cache
end

function Arrays.evaluate!(
  cache, k::InversePushforward, v_phys::Number, args...
)
  v_ref_basis, pf_cache = cache
  change = evaluate!(pf_cache,k.pushforward,v_ref_basis,args...)
  return v_phys⋅inv(change)
end

function evaluate!(
  cache, k::InversePushforward, f_phys::AbstractVector{<:Field}, args...
)
  Broadcasting(Operation(k))(f_phys,args...)
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
  ̂σ = F**(σ) := σ∘F* : V̂ -> R
"""
struct Pullback{PF} <: Map
  pushforward::PF
  function Pullback(pushforward::Pushforward)
    PF = typeof(pushforward)
    new{PF}(pushforward)
  end
end

function Arrays.lazy_map(
  ::typeof(evaluate),k::LazyArray{<:Fill{<:Pullback}},ref_cell_fields::AbstractArray
)
  pb = k.maps.value
  phys_cell_dofs, pf_args = k.args
  phys_cell_fields = lazy_map(pb.pushforward,ref_cell_fields,pf_args...)
  return lazy_map(evaluate,phys_cell_dofs,phys_cell_fields)
end

function evaluate!(
  cache, k::PullBack, σ_phys::MomentBasedDofBasis, args...
)
  Broadcasting(Operation(k))(f_phys,args...)
end

# InversePullback

"""
    struct InversePullback{PF <: Pushforward} <: Map end

Represents the inverse of the pullback map F**, defined as 
  (F**)^-1 : V̂* -> V*
where 
  - V̂* is a dof space on the reference cell K̂ and 
  - V* is a dof space on the physical cell K.
Its action on reference dofs ̂σ : V̂ -> R is defined in terms of the pushforward map F* as
  σ = F**(̂σ) := ̂σ∘(F*)^-1 : V -> R
"""
struct InversePullback{PF} <: Map
  pushforward::PF
  function InversePullback(pushforward::Pushforward)
    PF = typeof(pushforward)
    new{PF}(pushforward)
  end
end

Arrays.inverse_map(pb::Pullback) = InversePullback(pb.pushforward)
Arrays.inverse_map(ipb::InversePullback) = Pullback(ipb.pushforward)

function Arrays.lazy_map(
  ::typeof(evaluate),k::LazyArray{<:Fill{<:InversePullback}},phys_cell_fields::AbstractArray
)
  pb = inverse_map(k.maps.value)
  ref_cell_dofs, pf_args = k.args
  ref_cell_fields = lazy_map(inverse_map(pb.pushforward),phys_cell_fields,pf_args...)
  return lazy_map(evaluate,ref_cell_dofs,ref_cell_fields)
end

# ContraVariantPiolaMap

struct ContraVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::ContraVariantPiolaMap, v_ref::Number, Jt::Number, sign_flip::Bool
)
  sign  = (-1)^sign_flip
  idetJ = 1. / meas(Jt)
  return (sign*v_ref)⋅(idetJ*Jt)
end

# CoVariantPiolaMap

struct CoVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::CoVariantPiolaMap, v_ref::Number, Jt::Number
)
  # we right-multiply to compute the gradient correctly
  return v_ref⋅transpose(inv(Jt))
end
