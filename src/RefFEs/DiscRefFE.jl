
"""
Reference Finite Element without change-of-basis, all the DOFs belong to the
polytope itself
"""
struct DiscRefFE{D,T} <: RefFE{D,T}
  polytope::Polytope{D}
  shfbasis::Basis{D,T}
  nfacedofs::Vector{Vector{Int}}
  # this type is unstable
end

"""
Constructor that given a vector of orders, generates high-order Lagrangian
RefFEs on n-cubes or n-tets.
"""
function DiscRefFE(::Type{T}, p::Polytope{D}, orders::Vector{Int}) where {D,T}
  @assert length(orders) == D
  @assert D > 0
  # @santiagobadia : This part is missing for anisotropic orders !
  basis = Gridap.RefFEs._monomial_basis(p,T,orders[1])
  # aux = zeros(Float64,numlocaldofs(dofsb),numlocaldofs(dofsb))
  # @assert numlocaldofs(dofsb) == length(basis)
  return DiscRefFE{D,T}(p, basis, nfacedofs)
end

"""
Version of the constructor for a scalar order
"""
function DiscRefFE(::Type{T}, p::Polytope{D}, order::Int) where {D,T}
  _order = order*ones(Int,D)
  return DiscRefFE(T,p,_order)
end
