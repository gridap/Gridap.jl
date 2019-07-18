module RefFEs

using Gridap
using Gridap.Helpers
using Gridap.DOFBases: LagrangianDOFBasis

export RefFE
export LagrangianRefFE
export shfbasis
export polytope

import Gridap: evaluate, evaluate!

# Abstract types and interfaces

"""
Abstract Reference Finite Element
"""
abstract type RefFE{D,T} end

dofs(this::RefFE{D,T} where {D,T})::DOFBasis{D,T} = @abstractmethod

# permutation(this::RefFE, nf::Int, cell_vertex_gids::AbstractVector{Int},
# nface_vertex_gids::AbstractVector{Int}, nface_order::Int)::Vector{Int}
# = @abstractmethod
# @santiagobadia : To do in the future, not needed for the moment

polytope(this::RefFE{D,T} where {D,T})::Polytope{D} = @abstractmethod

shfbasis(this::RefFE{D,T} where {D,T})::Basis{D,T} = @abstractmethod

nfacedofs(this::RefFE{D,T} where {D,T})::Vector{Vector{Int}} = @abstractmethod

"""
Reference Finite Element a la Ciarlet, i.e., it relies on a local function
(polynomial) space, an array of nodes (DOFs), and a polytope (cell topology).
The rank of the approximating field can be arbitrary. The current implementation
relies on the prebasis (e.g., monomial basis of polynomials) and a
change-of-basis (using the node array) to generate the canonical basis, i.e.,
the shape functions.
"""
struct LagrangianRefFE{D,T} <: RefFE{D,T}
  polytope::Polytope{D}
  dofbasis::LagrangianDOFBasis{D,T}
  shfbasis::Basis{D,T}
  nfacedofs::Vector{Vector{Int}}
  # this type is unstable
end

function LagrangianRefFE{D,T}(polytope::Polytope{D},
  orders::Vector{Int64}) where {D,T}
  nodes=NodesArray(polytope,orders)
  dofsb = LagrangianDOFBasis{D,T}(nodes.coordinates)
  prebasis = MonomialBasis(T,orders)
  aux = zeros(Float64,numlocaldofs(dofsb),numlocaldofs(dofsb))
  @assert numlocaldofs(dofsb) == length(prebasis)
  changeofbasis=inv(evaluate!(dofsb,prebasis,aux))
  basis = change_basis(prebasis, changeofbasis)
  nfacedofs=nodes.nfacenodes
  LagrangianRefFE{D,T}(polytope, dofsb, basis, nfacedofs)
end

function _high_order_lagrangian_reffe(p::Polytope{D}, T, order) where {D}
  nodes, nfacedofs = _high_order_lagrangian_nodes_polytope(p,order)
  dofsb = Gridap.RefFEs.LagrangianDOFBasis{D,T}(nodes)
  prebasis = Gridap.RefFEs._monomial_basis(p,T,order)
  aux = zeros(Float64,numlocaldofs(dofsb),numlocaldofs(dofsb))
  @assert numlocaldofs(dofsb) == length(prebasis)
  changeofbasis=inv(evaluate!(dofsb,prebasis,aux))
  basis = change_basis(prebasis, changeofbasis)
  return LagrangianRefFE{D,T}(p, dofsb, basis, nfacedofs)
end

function _linear_lagrangian_reffe(polytope::Polytope{D},T) where {D}
  function _linear_nfacedofs(p)
    nfacedofs = Vector{Int}[]
    for (inf,nf) in enumerate(p.nfaces[p.nf_dim[end][1]])
      push!(nfacedofs, [inf])
    end
    for id in 2:length(p.nf_dim[end])
      for (inf,nf) in enumerate(p.nfaces[p.nf_dim[end][id]])
        push!(nfacedofs, Int[])
      end
    end
    return nfacedofs
  end
  nodes = Gridap.Polytopes.vertices_coordinates(polytope)
  dofsb = LagrangianDOFBasis{D,T}(nodes)
  order = 1
  prebasis = _monomial_basis(polytope,T,order)
  aux = zeros(Float64,numlocaldofs(dofsb),numlocaldofs(dofsb))
  @assert numlocaldofs(dofsb) == length(prebasis)
  changeofbasis=inv(evaluate!(dofsb,prebasis,aux))
  basis = change_basis(prebasis, changeofbasis)
  nfacedofs=_linear_nfacedofs(polytope)
  return LagrangianRefFE{D,T}(polytope, dofsb, basis, nfacedofs)
end

function _high_order_lagrangian_nodes_polytope(p::Polytope, order)
  vs_p = Gridap.Polytopes.vertices_coordinates(p)
  ns_float_p = [i.array for i in vs_p]
  ref_ps = Gridap.Polytopes._nface_ref_polytopes(p)
  rfe_p = Gridap.RefFEs._linear_lagrangian_reffe(p,Float64)
  nfacedofs = copy(rfe_p.nfacedofs)
  k = length(vs_p)
  for nf_dim = 1:length(p.nf_dim[end])-1
    nfs_dim = p.nf_dim[end][nf_dim+1]
    nfs_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,nf_dim,0)
    for (i_nf_dim,i_nf) in enumerate(p.nf_dim[end][nf_dim+1])
      ref_p = ref_ps[i_nf]
      _order = Tuple(order*ones(Int,length(ref_p.extrusion)))
      ns = Gridap.Polytopes.generate_interior_nodes(ref_p, _order)
      ns_float = Gridap.Polytopes._equidistant_nodes_coordinates(ns,_order)
      rfe = Gridap.RefFEs._linear_lagrangian_reffe(ref_p,Float64)
      nf_vs = nfs_vs[i_nf_dim]
      vs = vs_p[nf_vs]
      if ( length(ns_float) > 0 )
        a = evaluate(rfe.shfbasis,ns_float)
        nfacedofs[i_nf] = [k+i for i in 1:length(ns)]
        k += length(ns)
        hons = Gridap.RefFEs._map_high_order_lagrangian_nodes(a, vs)
        for i in hons
          push!(ns_float_p, i)
        end
      else
        nfacedofs[i_nf] = Int[]
      end
    end
  end
  return ns_float_p, nfacedofs
end

function _monomial_basis(p::Polytope{D}, T, order) where D
  if (_is_hex(p))
    orders = order*ones(Int,D)
    prebasis = MonomialBasis(T,orders)
  elseif (_is_tet(p))
    filter(e,order) = sum(e) <= order
    prebasis = MonomialBasis(D, T, filter, order)
  end
end

function _is_hex(p)
  for i in 2:length(p.extrusion)
    if p.extrusion[i] != 1 return false && break end
  end
  return true
end

function _is_tet(p)
  for i in 2:length(p.extrusion)
    if p.extrusion[i] != 2 return false && break end
  end
  return true
end

dofs(this::LagrangianRefFE{D,T} where {D,T}) = this.dofbasis

polytope(this::LagrangianRefFE{D,T} where {D,T}) = this.polytope

shfbasis(this::LagrangianRefFE{D,T} where {D,T}) = this.shfbasis

nfacedofs(this::LagrangianRefFE{D,T} where {D,T}) = this.nfacedofs

function _map_high_order_lagrangian_nodes(shfs_ns, vs)
  a = shfs_ns
  T = typeof(vs[1].array)
  ndofs, npoints = size(a)
  v = Vector{T}(undef,npoints)
  for j in 1:npoints
    aux = zero(T)
    for i in 1:ndofs
      aux += outer(a[i,j],vs[i]).array
    end
    v[j] = aux
  end
  return v
end

end # module RefFEs
