
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
  ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{<:Pushforward}}}}
)
  cell_ref_basis, args = a.args
  cell_ref_gradient = lazy_map(Broadcasting(∇),cell_ref_basis)
  return lazy_map(a.maps.value,cell_ref_gradient,args...)
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

# Pushforward

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
  ::typeof{evaluate},k::LazyArray{<:Fill{<:Pushforward}},ref_cell_basis
)
  pb = k.maps.value
  phys_cell_dofs, cell_map, pf_args = k.args
  phys_cell_basis = lazy_map(pb.pushforward,ref_cell_basis,cell_map,pf_args...)
  return lazy_map(evaluate,phys_cell_dofs,phys_cell_basis)
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
  ::typeof{evaluate},k::LazyArray{<:Fill{<:InversePullback}},phys_cell_basis
)
  pb = inverse_map(k.maps.value)
  ref_cell_dofs, cell_map, pf_args = k.args
  ref_cell_basis = lazy_map(inverse_map(pb.pushforward),phys_cell_basis,cell_map,pf_args...)
  return lazy_map(evaluate,ref_cell_dofs,ref_cell_basis)
end

# ContraVariantPiolaMap

struct ContraVariantPiolaMap <: Pushforward end

function lazy_map(
  k::ContraVariantPiolaMap,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
  cell_map::AbstractArray{<:Field},
  sign_flip::AbstractArray{<:AbstractArray{<:Field}}
)
  cell_Jt = lazy_map(∇,cell_map)
  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_Jt,sign_flip)
end

function evaluate!(
  cache,::ContraVariantPiolaMap,
  v::Number,
  Jt::Number,
  sign_flip::Bool
)
  idetJ = 1/meas(Jt)
  ((-1)^sign_flip*v)⋅(idetJ*Jt)
end
