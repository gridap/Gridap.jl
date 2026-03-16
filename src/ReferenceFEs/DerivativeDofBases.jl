struct PointDerivativeValue{P} <: Dof
  point::P
  direction::P
end

"""
    struct PointDerivativeDofBasis{P,V} <: AbstractArray{PointDerivativeValue{P}}
      nodes::Vector{P}
      dof_to_node::Vector{Int}
      dof_to_comp::Vector{Int}
      node_and_comp_to_dof::Vector{V}
    end

Type that implements a PointDerivative DOF basis, that is a [`Dof`](@ref) basis
where the DOF are nodal evaluations of directional derivatives. If the value
type `V` is not scalar ([`::MultiValue`](@ref MultiValue)'d), several DoFs may
be associated to the same node, one for each independant component of `V`.

Fields:

- `nodes` vector of points ([`P<:Point`](@ref Point)) storing the nodal coordinates
- `node_and_comp_to_dof` vector such that `node_and_comp_to_dof[node][comp]` returns the dof associated with node `node` and the component `comp` in the type `V`.
- `dof_to_node` vector of integers such that `dof_to_node[dof]` returns the node id associated with dof id `dof`.
- `dof_to_comp` vector of integers such that `dof_to_comp[dof]` returns the component id associated with dof id `dof`.
"""
struct PointDerivativeDofBasis{P,V} <: AbstractVector{PointDerivativeValue{P}}
  nodes::Vector{P}
  dof_to_node::Vector{Int}
  dof_to_comp::Vector{Int}
  node_and_comp_to_dof::Vector{V}
end

Base.size(a::PointDerivativeDofBasis) = (length(a.dof_to_node),)
Base.axes(a::PointDerivativeDofBasis) = (axes(a.dof_to_node,1),)
# @santiagobadia : Not sure we want to create the monomial machinery
Base.getindex(a::PointDerivativeDofBasis,i::Integer) = PointDerivativeValue(a.nodes[a.dof_to_node[i]])
Base.IndexStyle(::PointDerivativeDofBasis) = IndexLinear()

# This one takes a basis and replaces the nodes
function PointDerivativeDofBasis(dofs::PointDerivativeDofBasis{P},nodes::Vector{P}) where P
  @check length(nodes) == length(dofs.nodes)
  PointDerivativeDofBasis(
    nodes,
    dofs.dof_to_node,
    dofs.dof_to_comp,
    dofs.node_and_comp_to_dof)
end

"""
    PointDerivativeDofBasis(::Type{V}, nodes::Vector{<:Point})

Creates a `PointDerivativeDofBasis` for fields of value type `V` associated
with the vector of nodal coordinates `nodes`.
"""
function PointDerivativeDofBasis(::Type{V},nodes::Vector{<:Point}) where V
  r = _generate_dof_layout_node_major(V,length(nodes))
  PointDerivativeDofBasis(nodes,r...)
end

"""
    get_nodes(b::PointDerivativeDofBasis)

Get the vector of DoF nodes of `b`.
"""
get_nodes(b::PointDerivativeDofBasis) = b.nodes
get_dof_to_node(b::PointDerivativeDofBasis) = b.dof_to_node
get_dof_to_comp(b::PointDerivativeDofBasis) = b.dof_to_comp

#function _generate_dof_layout_node_major(::Type{<:Real},nnodes::Integer)
#  ndofs = nnodes
#  dof_to_comp = ones(Int,ndofs)
#  dof_to_node = collect(1:nnodes)
#  node_and_comp_to_dof = collect(1:ndofs)
#  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
#end
#
## Node major implementation
#function _generate_dof_layout_node_major(::Type{V},nnodes::Integer) where V<:MultiValue
#  Vi = change_eltype(V,Int)
#  ncomps = num_indep_components(Vi)
#  ndofs = ncomps*nnodes
#  dof_to_comp = zeros(Int,ndofs)
#  dof_to_node = zeros(Int,ndofs)
#  node_and_comp_to_dof = Vector{Vi}(undef,nnodes)
#  m = zero(MVector{ncomps,Int})
#  for node in 1:nnodes
#    for comp in 1:ncomps
#      o = nnodes*(comp-1)
#      dof = node+o
#      dof_to_comp[dof] = comp
#      dof_to_node[dof] = node
#      m[comp] = dof
#    end
#    node_and_comp_to_dof[node] = Tuple(m)
#  end
#  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
#end

function return_cache(b::PointDerivativeDofBasis,field)
  cf = return_cache(field,b.nodes)
  vals = evaluate!(cf,field,b.nodes)
  ndofs = length(b.dof_to_node)
  r = _lagr_dof_cache(vals,ndofs)
  c = CachedArray(r)
  (c, cf)
end

#function _lagr_dof_cache(node_comp_to_val::AbstractVector,ndofs)
#  V = eltype(node_comp_to_val)
#  r = zeros(eltype(V),ndofs)
#end
#
#function _lagr_dof_cache(node_pdof_comp_to_val::AbstractMatrix,ndofs)
#  _, npdofs = size(node_pdof_comp_to_val)
#  V = eltype(node_pdof_comp_to_val)
#  r = zeros(eltype(V),ndofs,npdofs)
#end

function evaluate!(cache,b::PointDerivativeDofBasis,field)
  c, cf = cache
  vals = evaluate!(cf,field,b.nodes)
  ndofs = length(b.dof_to_node)
  V = eltype(vals)
  ncomps = num_indep_components(V)
  @check ncomps == num_indep_components(eltype(b.node_and_comp_to_dof)) """\n
  Unable to evaluate PointDerivativeDofBasis. The number of components of the
  given Field does not match with the PointDerivativeDofBasis.

  If you are trying to interpolate a function on a FESpace make sure that
  both objects have the same value type.

  For instance, trying to interpolate a vector-valued function on a scalar-valued FE space
  would raise this error.
  """
  _evaluate_lagr_dof!(c,vals,b.node_and_comp_to_dof,ndofs,ncomps)
end


#function _evaluate_lagr_dof!(c::AbstractVector,node_comp_to_val,node_and_comp_to_dof,ndofs,ncomps)
#  setsize!(c,(ndofs,))
#  r = c.array
#  for node in LinearIndices(node_and_comp_to_dof)
#    comp_to_dof = node_and_comp_to_dof[node]
#    comp_to_val = node_comp_to_val[node]
#    for comp in 1:ncomps
#      dof = indep_comp_getindex(comp_to_dof,comp)
#      val = indep_comp_getindex(comp_to_val,comp)
#      r[dof] = val
#    end
#  end
#  r
#end
#
#function _evaluate_lagr_dof!(c::AbstractMatrix,node_pdof_comp_to_val,node_and_comp_to_dof,ndofs,ncomps)
#  _, npdofs = size(node_pdof_comp_to_val)
#  setsize!(c,(ndofs,npdofs))
#  r = c.array
#  for node in LinearIndices(node_and_comp_to_dof)
#    comp_to_dof = node_and_comp_to_dof[node]
#    for pdof in 1:npdofs
#      comp_to_val = node_pdof_comp_to_val[node,pdof]
#      for comp in 1:ncomps
#        dof = indep_comp_getindex(comp_to_dof,comp)
#        val = indep_comp_getindex(comp_to_val,comp)
#        r[dof,pdof] = val
#      end
#    end
#  end
#  r
#end
