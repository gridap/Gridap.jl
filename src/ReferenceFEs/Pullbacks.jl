
"""
    abstract type Pullback <: Map end

Represents a pullback map F*, defined as 
  F* : V̂ -> V
where 
  - V̂ is a function space on the reference cell K̂ and 
  - V is a function space on the physical cell K.
"""
abstract type Pullback <: Map end

function lazy_map(
  ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{<:Pullback}}}}
)
  cell_ref_basis, args = a.args
  cell_ref_gradient = lazy_map(Broadcasting(∇),cell_ref_basis)
  return lazy_map(a.maps.value,cell_ref_gradient,args...)
end

# InversePullback

"""
    struct InversePullback{PB <: Pullback} <: Map end

Represents the inverse of a pullback map F*, defined as 
  (F*)^-1 : V -> V̂
where 
  - V̂ is a function space on the reference cell K̂ and 
  - V is a function space on the physical cell K.
"""
struct InversePullback{PB} <: Map
  pullback::PB
  function InversePullback(pullback::Pullback)
    PB = typeof(pullback)
    new{PB}(pullback)
  end
end

Arrays.inverse_map(pb::Pullback) = InversePullback(pb)
Arrays.inverse_map(ipb::InversePullback) = ipb.pullback

# Pushforward

"""
    struct Pushforward{PB <: Pullback} <: Map end

Represents a pushforward map F**, defined as 
  F** : V* -> V̂*
where 
  - V̂* is a dof space on the reference cell K̂ and 
  - V* is a dof space on the physical cell K.
Its action on physical dofs σ : V -> R is defined in terms of the pullback map F* as
  F**(σ) := σ∘F* : V̂ -> R
"""
struct Pushforward{PB} <: Map
  pullback::PB
  function Pushforward(pullback::Pullback)
    PB = typeof(pullback)
    new{PB}(pullback)
  end
end

function Arrays.lazy_map(
  ::typeof{evaluate},k::LazyArray{<:Fill{<:Pushforward}},ref_cell_basis
)
  pf = k.maps.value
  phys_cell_dofs, cell_map, pb_args = k.args
  phys_cell_basis = lazy_map(pf.pullback,ref_cell_basis,cell_map,pb_args...)
  return lazy_map(evaluate,phys_cell_dofs,phys_cell_basis)
end

# InversePushforward

"""
    struct InversePushforward{PB <: Pullback} <: Map end

Represents the inverse of the pushforward map F**, defined as 
  (F**)^-1 : V̂* -> V*
where 
  - V̂* is a dof space on the reference cell K̂ and 
  - V* is a dof space on the physical cell K.
Its action on reference dofs ̂σ : V -> R is defined in terms of the pullback map F* as
  F**(̂σ) := ̂σ∘(F*)^-1 : V -> R
"""
struct InversePushforward{PB} <: Map
  pullback::PB
  function InversePushforward(pullback::Pullback)
    PB = typeof(pullback)
    new{PB}(pullback)
  end
end

Arrays.inverse_map(pf::Pushforward) = InversePushforward(pf.pullback)
Arrays.inverse_map(ipf::InversePushforward) = Pushforward(ipf.pullback)

function Arrays.lazy_map(
  ::typeof{evaluate},k::LazyArray{<:Fill{<:InversePushforward}},phys_cell_basis
)
  pf = inverse_map(k.maps.value)
  ref_cell_dofs, cell_map, pb_args = k.args
  ref_cell_basis = lazy_map(inverse_map(pf.pullback),phys_cell_basis,cell_map,pb_args...)
  return lazy_map(evaluate,ref_cell_dofs,ref_cell_basis)
end

# ContraVariantPiolaMap

struct ContraVariantPiolaMap <: Pullback end

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
