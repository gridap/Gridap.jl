"""
Thoughts on this implementation: 

- In order to be able to dispatch on the choice of integration space, we could 
  use a Val{Symbol} as a emplate parameter. For instance 

    struct Myobj{A}
      function Myobj(sym::Symbol)
        A = typeof(Val(sym))
        new{A}()
      end
    end

    and then dispatch mul!() depending on the symbol (:from/:to).

- However, this might be limiting / not general enough. For instance, when we do 
  not have uniform refinement (some cells get refined, other get coarsened) we will
  have an integration space which is a combination of 'from' and 'to'. 
  We could solve this with a third symbol dispatch like :mixed, but it is not very elegant. 

- Another solution would be to always apply 'change_domain' to both the FEFunction and the 
  test FEBasis. This is general since all changes are done in the background. This might be 
  the most elegant option, but is it the most efficient? 
"""

"""
TODOs: 
  - Implement change_domain(FEFunction,RefinedTriangulation)
  - Implement change_domain(FEBasis,RefinedTriangulation)
  - Implement a choice for the projection
"""

struct RefinementTransferOperator{T,A,B,C} <: AbstractMatrix{T}
  from   ::A
  to     ::B
  caches ::C

  function RefinementTransferOperator(T,from,to,caches)
    A = typeof(from)
    B = typeof(to)
    C = typeof(caches)
    new{T,A,B,C}(from,to,caches)
  end
end

function RefinementTransferOperator(from::FESpace,to::FESpace; qdegree=1)
  @assert isa(from,TrialFESpace)
  @assert isa(to,TrialFESpace)

  Ω_from = get_triangulation(from)
  Ω_to   = get_triangulation(to)
  @assert isa(Ω_from,RefinedTriangulation) || isa(Ω_to,RefinedTriangulation)
  @assert is_change_possible(Ω_from,Ω_to)

  # Choose integration space (finest)
  Ω  = best_target(Ω_from,Ω_to)
  dΩ = Measure(Ω,qdegree)
  U  = (Ω === Ω_from) ? from : to
  V  = U.space
  vh = change_domain(get_fe_basis(V),Ω)

  # Prepare system. TODO: Choosing the projection method should be left to the user. 
  sysmat = assemble_mass_matrix(Ω_to,to,to.space,qdegree)
  sysvec = zeros(size(sysmat,1))
  assem  = SparseMatrixAssembler(to,from.space)
  rhs(uh,vh) = ∫(vh⋅uh) * dΩ

  cache = sysmat, sysvec, rhs, assem, Ω, dΩ, U, V, vh
  return RefinementTransferOperator(eltype(sysmat),from,to,cache)
end

# Solves the problem (uh,vh)_to = (uh_from,vh)_Ω for all vh in Vh_to
function mul!(y,A::RefinementTransferOperator,x)
  sysmat, sysvec, rhs, assem, Ω, dΩ, U, V , vh_Ω = A.cache

  # Bring uh to the integration domain
  uh_from = FEFunction(A.from,x)
  uh_Ω    = change_domain(uh_from,Ω)

  # Assemble rhs vector
  contr   = rhs(vh_Ω,uh_Ω)
  vecdata = collect_cell_vector(A.to.space,contr)
  assemble_vector!(sysvec,assem,vecdata)

  # Solve projection
  cg!(y,sysmat,sysvec)

  return y
end

function assemble_mass_matrix(Ω,Uh,Vh,qdegree)
  dΩ = Measure(Ω,qdegree)
  a(u,v) = ∫(v⋅u)*dΩ
  return assemble_matrix(a,Uh,Vh)
end