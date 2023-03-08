
"""
"""
struct FineToCoarseDofBasis{T,A,B,C} <: AbstractVector{T}
  dof_basis :: A
  rrule     :: B
  child_ids :: C

  function FineToCoarseDofBasis(dof_basis::AbstractVector{T},rrule::RefinementRule) where {T<:Dof}
    nodes = get_nodes(dof_basis)
    child_ids = map(x -> x_to_cell(rrule,x),nodes)

    A = typeof(dof_basis)
    B = typeof(rrule)
    C = typeof(child_ids)
    new{T,A,B,C}(dof_basis,rrule,child_ids)
  end
end

Base.size(a::FineToCoarseDofBasis) = size(a.dof_basis)
Base.axes(a::FineToCoarseDofBasis) = axes(a.dof_basis)
Base.getindex(a::FineToCoarseDofBasis,i::Integer) = getindex(a.dof_basis,i)
Base.IndexStyle(a::FineToCoarseDofBasis) = IndexStyle(a.dof_basis)

ReferenceFEs.get_nodes(a::FineToCoarseDofBasis) = get_nodes(a.dof_basis)

# Default behaviour
Arrays.return_cache(b::FineToCoarseDofBasis,field) = return_cache(b.dof_basis,field)
Arrays.evaluate!(cache,b::FineToCoarseDofBasis,field) = evaluate!(cache,b.dof_basis,field)

# Spetialized behaviour
function Arrays.return_cache(s::FineToCoarseDofBasis{T,<:LagrangianDofBasis},field::FineToCoarseField) where T
  b     = s.dof_basis
  cf    = return_cache(field,b.nodes,s.child_ids)
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  ndofs = length(b.dof_to_node)
  r = ReferenceFEs._lagr_dof_cache(vals,ndofs)
  c = CachedArray(r)
  return (c, cf)
end

function Arrays.evaluate!(cache,s::FineToCoarseDofBasis{T,<:LagrangianDofBasis},field::FineToCoarseField) where T
  c, cf = cache
  b     = s.dof_basis
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  ndofs = length(b.dof_to_node)
  T2    = eltype(vals)
  ncomps = num_components(T2)
  @check ncomps == num_components(eltype(b.node_and_comp_to_dof)) """\n
  Unable to evaluate LagrangianDofBasis. The number of components of the
  given Field does not match with the LagrangianDofBasis.

  If you are trying to interpolate a function on a FESpace make sure that
  both objects have the same value type.

  For instance, trying to interpolate a vector-valued function on a scalar-valued FE space
  would raise this error.
  """
  ReferenceFEs._evaluate_lagr_dof!(c,vals,b.node_and_comp_to_dof,ndofs,ncomps)
end

function Arrays.return_cache(s::FineToCoarseDofBasis{T,<:MomentBasedDofBasis},field::FineToCoarseField) where T
  b     = s.dof_basis
  cf    = return_cache(field,b.nodes,s.child_ids)
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  ndofs = num_dofs(b)
  r = ReferenceFEs._moment_dof_basis_cache(vals,ndofs)
  c = CachedArray(r)
  return (c, cf)
end

function Arrays.evaluate!(cache,s::FineToCoarseDofBasis{T,<:MomentBasedDofBasis},field::FineToCoarseField) where T
  c, cf = cache
  b     = s.dof_basis
  vals  = evaluate!(cf,field,b.nodes,s.child_ids)
  dofs  = c.array
  ReferenceFEs._eval_moment_dof_basis!(dofs,vals,b)
  dofs
end


"""
  Wrapper for a ReferenceFE which is specialised for 
  efficiently evaluating FineToCoarseFields. 
"""
struct FineToCoarseRefFE{T,D,A} <: ReferenceFE{D}
  reffe     :: T
  dof_basis :: A

  function FineToCoarseRefFE(reffe::ReferenceFE{D},dof_basis::FineToCoarseDofBasis) where D
    T = typeof(reffe)
    A = typeof(dof_basis)
    new{T,D,A}(reffe,dof_basis)
  end
end

ReferenceFEs.num_dofs(reffe::FineToCoarseRefFE)      = num_dofs(reffe.reffe)
ReferenceFEs.get_polytope(reffe::FineToCoarseRefFE)  = get_polytope(reffe.reffe)
ReferenceFEs.get_prebasis(reffe::FineToCoarseRefFE)  = get_prebasis(reffe.reffe)
ReferenceFEs.get_dof_basis(reffe::FineToCoarseRefFE) = reffe.dof_basis
ReferenceFEs.Conformity(reffe::FineToCoarseRefFE)    = Conformity(reffe.reffe)
ReferenceFEs.get_face_dofs(reffe::FineToCoarseRefFE) = get_face_dofs(reffe.reffe)
ReferenceFEs.get_shapefuns(reffe::FineToCoarseRefFE) = get_shapefuns(reffe.reffe)
ReferenceFEs.get_metadata(reffe::FineToCoarseRefFE)  = get_metadata(reffe.reffe)
ReferenceFEs.get_orders(reffe::FineToCoarseRefFE)    = get_orders(reffe.reffe)
ReferenceFEs.get_order(reffe::FineToCoarseRefFE)     = get_order(reffe.reffe)

ReferenceFEs.Conformity(reffe::FineToCoarseRefFE,sym::Symbol) = Conformity(reffe.reffe,sym)
ReferenceFEs.get_face_own_dofs(reffe::FineToCoarseRefFE,conf::Conformity) = get_face_own_dofs(reffe.reffe,conf)


function ReferenceFEs.ReferenceFE(p::Polytope,rrule::RefinementRule,name::ReferenceFEName,order)
  FineToCoarseRefFE(p,rrule,name,Float64,order)
end

function ReferenceFEs.ReferenceFE(p::Polytope,rrule::RefinementRule,name::ReferenceFEName,::Type{T},order) where T
  FineToCoarseRefFE(p,rrule,name,T,order)
end

function FineToCoarseRefFE(p::Polytope,rrule::RefinementRule,name::ReferenceFEName,::Type{T},order) where T
  @check p == get_polytope(rrule)
  reffe = ReferenceFE(p,name,T,order)
  dof_basis = FineToCoarseDofBasis(get_dof_basis(reffe),rrule)
  return FineToCoarseRefFE(reffe,dof_basis)
end

# FESpaces constructors

function FESpaces.TestFESpace(model::DiscreteModel,rrules::AbstractVector{<:RefinementRule},reffe::Tuple{<:ReferenceFEName,Any,Any};kwargs...)
  @check num_cells(model) == length(rrules)
  @check all(CompressedArray(get_polytopes(model),get_cell_type(model)) .== lazy_map(get_polytope,rrules))
  basis, reffe_args, reffe_kwargs = reffe
  reffes = lazy_map(rr -> ReferenceFE(get_polytope(rr),rr,basis,reffe_args...;reffe_kwargs...),rrules)
  return TestFESpace(model,reffes;kwargs...)
end