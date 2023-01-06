
"""
"""
struct FineToCoarseDofBasis{A<:AbstractVector{Dof},B<:RefinementRule,C} <: AbstractVector{Dof}
  dof_basis :: A
  rrule     :: B
  child_ids :: C

  function FineToCoarseDofBasis(dof_basis::AbstractVector{Dof},rrule::RefinementRule)
    child_ids = map(rrule.x_to_cell,dof_basis.nodes)

    A = typeof(dof_basis)
    B = typeof(rrule)
    C = typeof(child_ids)
    return new{A,B,C}(dof_basis,rrule,child_ids)
  end
end

# Default behaviour
get_cache(b::FineToCoarseDofBasis,field) = return_cache(b.dof_basis,field)
evaluate!(cache,b::FineToCoarseDofBasis,field) = evaluate!(cache,b.dof_basis,field)

# Spetialised behaviour
function get_cache(s::FineToCoarseDofBasis{<:LagrangianDofBasis},field::FineToCoarseField)
  b     = s.dof_basis
  cf    = return_cache(field,b.nodes,s.child_ids)
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  ndofs = length(b.dof_to_node)
  r = _lagr_dof_cache(vals,ndofs)
  c = CachedArray(r)
  return (c, cf)
end

function evaluate!(cache,s::FineToCoarseDofBasis{<:LagrangianDofBasis},field::FineToCoarseField)
  c, cf = cache
  b     = s.dof_basis
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  ndofs = length(b.dof_to_node)
  T = eltype(vals)
  ncomps = num_components(T)
  @check ncomps == num_components(eltype(b.node_and_comp_to_dof)) """\n
  Unable to evaluate LagrangianDofBasis. The number of components of the
  given Field does not match with the LagrangianDofBasis.

  If you are trying to interpolate a function on a FESpace make sure that
  both objects have the same value type.

  For instance, trying to interpolate a vector-valued funciton on a scalar-valued FE space
  would raise this error.
  """
  _evaluate_lagr_dof!(c,vals,b.node_and_comp_to_dof,ndofs,ncomps)
end

function get_cache(s::FineToCoarseDofBasis{<:MomentBasedDofBasis},field::FineToCoarseField)
  b     = s.dof_basis
  cf    = return_cache(field,b.nodes,s.child_ids)
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  ndofs = num_dofs(b)
  r = _moment_dof_basis_cache(vals,ndofs)
  c = CachedArray(r)
  return (c, cf)
end

function evaluate!(cache,s::FineToCoarseDofBasis{<:MomentBasedDofBasis},field)
  c, cf = cache
  b     = s.dof_basis
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  dofs  = c.array
  _eval_moment_dof_basis!(dofs,vals,b)
  dofs
end


"""
  Wrapper for a ReferenceFE which is spetialised for 
  efficiently evaluating FineToCoarseFields. 
"""
struct FineToCoarseRefFE{T,D,A} <: ReferenceFE{D}
  reffe     :: T
  dof_basis :: A

  function FineToCoarseRefFE(reffe::ReferenceFE{D},dof_basis::FineToCoarseDofBasis) where D
    T = typeof(reffe)
    A = typeof(dof_basis)
    return new{T,A}(reffe,dof_basis)
  end
end

num_dofs(reffe::FineToCoarseRefFE)      = num_dofs(reffe.reffe)
get_polytope(reffe::FineToCoarseRefFE)  = get_polytope(reffe.reffe)
get_prebasis(reffe::FineToCoarseRefFE)  = get_prebasis(reffe.reffe)
get_dof_basis(reffe::FineToCoarseRefFE) = reffe.dof_basis
Conformity(reffe::FineToCoarseRefFE)    = Conformity(reffe.reffe)
get_face_dofs(reffe::FineToCoarseRefFE) = get_face_dofs(reffe.reffe)
get_shapefuns(reffe::FineToCoarseRefFE) = get_shapefuns(reffe.reffe)
get_metadata(reffe::FineToCoarseRefFE)  = get_metadata(reffe.reffe)

function ReferenceFE(p::Polytope,rrule::RefinementRule,name::ReferenceFEName,order)
  FineToCoarseRefFE(p,rrule,name,Float64,order)
end

function ReferenceFE(p::Polytope,rrule::RefinementRule,name::ReferenceFEName,::Type{T},order) where T
  FineToCoarseRefFE(p,rrule,name,T,order)
end

function FineToCoarseRefFE(p::Polytope,rrule::RefinementRule,name::ReferenceFEName,::Type{T},order) where T
  @check p == get_polytope(rrule)
  reffe = ReferenceFE(p,name,T,order)
  dof_basis = FineToCoarseDofBasis(get_dof_basis(reffe),rrule)
  return FineToCoarseRefFE(reffe,dof_basis)
end
