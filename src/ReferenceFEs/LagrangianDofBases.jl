
struct PointValue{P} <: Dof
  point::P
end

"""
    struct LagrangianDofBasis{P,V} <: AbstractArray{<:Dof}
      nodes::Vector{P}
      dof_to_node::Vector{Int}
      dof_to_comp::Vector{Int}
      node_and_comp_to_dof::Vector{V}
    end

Type that implements a Lagrangian dof basis.

Fields:

- `nodes::Vector{P}` vector of points (`P<:Point`) storing the nodal coordinates
- `node_and_comp_to_dof::Vector{V}` vector such that `node_and_comp_to_dof[node][comp]` returns the dof associated with node `node` and the component `comp` in the type `V`.
- `dof_to_node::Vector{Int}` vector of integers such that `dof_to_node[dof]` returns the node id associated with dof id `dof`.
- `dof_to_comp::Vector{Int}` vector of integers such that `dof_to_comp[dof]` returns the component id associated with dof id `dof`.

"""
struct LagrangianDofBasis{P,V} <: AbstractVector{PointValue{P}}
  nodes::Vector{P}
  dof_to_node::Vector{Int}
  dof_to_comp::Vector{Int}
  node_and_comp_to_dof::Vector{V}
end

@inline Base.size(a::LagrangianDofBasis) = (length(a.nodes),)
@inline Base.axes(a::LagrangianDofBasis) = (axes(a.nodes,1),)
# @santiagobadia : Not sure we want to create the monomial machinery
@inline Base.getindex(a::LagrangianDofBasis,i::Integer) = PointValue(a.nodes[i])
@inline Base.IndexStyle(::LagrangianDofBasis) = IndexLinear()

# This one takes a basis and replaces the nodes
function LagrangianDofBasis(dofs::LagrangianDofBasis{P},nodes::Vector{P}) where P
  @check length(nodes) == length(dofs.nodes)
  LagrangianDofBasis(
    nodes,
    dofs.dof_to_node,
    dofs.dof_to_comp,
    dofs.node_and_comp_to_dof)
end

"""
    LagrangianDofBasis(::Type{T},nodes::Vector{<:Point}) where T

Creates a `LagrangianDofBasis` for fields of value type `T` associated
with the vector of nodal coordinates `nodes`.
"""
function LagrangianDofBasis(::Type{T},nodes::Vector{<:Point}) where T
  r = _generate_dof_layout_node_major(T,length(nodes))
  LagrangianDofBasis(nodes,r...)
end

get_nodes(b::LagrangianDofBasis) = b.nodes

function _generate_dof_layout_node_major(::Type{<:Real},nnodes::Integer)
  ndofs = nnodes
  dof_to_comp = ones(Int,ndofs)
  dof_to_node = collect(1:nnodes)
  node_and_comp_to_dof = collect(1:ndofs)
  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

# Node major implementation
function _generate_dof_layout_node_major(::Type{T},nnodes::Integer) where T<:MultiValue
  ncomps = num_components(T)
  V = change_eltype(T,Int)
  ndofs = ncomps*nnodes
  dof_to_comp = zeros(Int,ndofs)
  dof_to_node = zeros(Int,ndofs)
  node_and_comp_to_dof = zeros(V,nnodes)
  m = zero(Mutable(V))
  for node in 1:nnodes
    for comp in 1:ncomps
      o = nnodes*(comp-1)
      dof = node+o
      dof_to_comp[dof] = comp
      dof_to_node[dof] = node
      m[comp] = dof
    end
    node_and_comp_to_dof[node] = m
  end
  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

function return_cache(b::LagrangianDofBasis,field)
  cf = return_cache(field,b.nodes)
  vals = evaluate!(cf,field,b.nodes)
  ndofs = length(b.dof_to_node)
  r = _lagr_dof_cache(vals,ndofs)
  c = CachedArray(r)
  (c, cf)
end

function _lagr_dof_cache(node_comp_to_val::AbstractVector,ndofs)
  T = eltype(node_comp_to_val)
  r = zeros(eltype(T),ndofs)
end

function _lagr_dof_cache(node_pdof_comp_to_val::AbstractMatrix,ndofs)
  _, npdofs = size(node_pdof_comp_to_val)
  T = eltype(node_pdof_comp_to_val)
  r = zeros(eltype(T),ndofs,npdofs)
end

@inline function evaluate!(cache,b::LagrangianDofBasis,field)
  c, cf = cache
  vals = evaluate!(cf,field,b.nodes)
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

function _evaluate_lagr_dof!(c::AbstractVector,node_comp_to_val,node_and_comp_to_dof,ndofs,ncomps)
  setsize!(c,(ndofs,))
  r = c.array
  for node in LinearIndices(node_and_comp_to_dof)
    comp_to_dof = node_and_comp_to_dof[node]
    comp_to_val = node_comp_to_val[node]
    for comp in 1:ncomps
      dof = comp_to_dof[comp]
      val = comp_to_val[comp]
      r[dof] = val
    end
  end
  r
end

function _evaluate_lagr_dof!(c::AbstractMatrix,node_pdof_comp_to_val,node_and_comp_to_dof,ndofs,ncomps)
  _, npdofs = size(node_pdof_comp_to_val)
  setsize!(c,(ndofs,npdofs))
  r = c.array
  for node in LinearIndices(node_and_comp_to_dof)
    comp_to_dof = node_and_comp_to_dof[node]
    for pdof in 1:npdofs
      comp_to_val = node_pdof_comp_to_val[node,pdof]
      for comp in 1:ncomps
        dof = comp_to_dof[comp]
        val = comp_to_val[comp]
        r[dof,pdof] = val
      end
    end
  end
  r
end
