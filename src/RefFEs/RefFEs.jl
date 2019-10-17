module RefFEs

using Gridap
using Gridap.Helpers
using Gridap.DOFBases: LagrangianDOFBasis

export RefFE
export LagrangianRefFE
export shfbasis
export polytope
export nfacedofs
export dofbasis
export field_type

import Gridap: evaluate, evaluate!

# Abstract types and interfaces

"""
Abstract Reference Finite Element
"""
abstract type RefFE{D,T} end

dofbasis(this::RefFE{D,T} where {D,T})::DOFBasis{D,T} = @abstractmethod

polytope(this::RefFE{D,T} where {D,T})::Polytope{D} = @abstractmethod

shfbasis(this::RefFE{D,T} where {D,T})::Basis{D,T} = @abstractmethod

nfacedofs(this::RefFE{D,T} where {D,T})::Vector{Vector{Int}} = @abstractmethod

function field_type(this::RefFE{D,T} ) where {D,T}
  return T
end

# @fverdugo this abstract interface needs a tester

"""
Extract the lists of dofs, but only for the nfaces of dimension dim
"""
function nfacedofs(reffe::RefFE{D},dim::Integer) where D
  @assert 0 <= dim
  @assert dim <= D
  nface_to_dofs = nfacedofs(reffe)
  p = polytope(reffe)
  nface_to_dofs[nf_dim(p)[end][dim+1]]
end

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
  nfacenodes::Vector{Vector{Int}}
  # this type is unstable
end

"""
Constructor that given a vector of orders, generates high-order Lagrangian
RefFEs on n-cubes or n-tets.
"""
function LagrangianRefFE(::Type{T}, p::Polytope{D}, orders::Vector{Int}) where {D,T}
  @assert length(orders) == D
  @assert D > 0
  if !(all(extrusion(p).array.== HEX_AXIS) || all(orders.==orders[1]))
    error("One can consider anisotropic orders on n-cubes only")
  end
  if (all(orders.==1))
    nodes, nfacenodes = _linear_lagrangian_nodes_polytope(p)
  else
    nodes, nfacenodes = _high_order_lagrangian_nodes_polytope(p,orders)
  end
  dofsb = Gridap.RefFEs.LagrangianDOFBasis(T,nodes)
  # @santiagobadia : This part is missing for anisotropic orders !
  prebasis = Gridap.RefFEs._monomial_basis(p,T,orders[1])
  aux = zeros(Float64,numlocaldofs(dofsb),numlocaldofs(dofsb))
  @assert numlocaldofs(dofsb) == length(prebasis)
  evaluate!(dofsb,prebasis,aux)
  changeofbasis=inv(aux)
  basis = change_basis(prebasis, changeofbasis)
  nfacedofs = _genereate_nface_to_dofs(nfacenodes,dofsb.node_and_comp_to_dof)
  return LagrangianRefFE{D,T}(p, dofsb, basis, nfacedofs, nfacenodes)
end

"""
Version of the constructor for a scalar order
"""
function LagrangianRefFE(::Type{T}, p::Polytope{D}, order::Int) where {D,T}
  _order = order*ones(Int,D)
  return LagrangianRefFE(T,p,_order)
end


dofbasis(this::LagrangianRefFE{D,T} where {D,T}) = this.dofbasis

polytope(this::LagrangianRefFE{D,T} where {D,T}) = this.polytope

shfbasis(this::LagrangianRefFE{D,T} where {D,T}) = this.shfbasis

nfacedofs(this::LagrangianRefFE{D,T} where {D,T}) = this.nfacedofs

function closurenfacedofs(this::LagrangianRefFE{D,T} where {D,T})
  cv1 = CellValueFromArray(this.polytope.nf_nfs) # cell to index
  cv2 = CellValueFromArray(this.nfacedofs)
  return collect(CellVectorByComposition(cv1,cv2))
end

# Generate the linear nodes by computing the polytope vertices. Only the
# vertices have nodes.
function _linear_lagrangian_nodes_polytope(p::Polytope)
  function _linear_nfacedofs(p)
    nfacedofs = Vector{Int}[]
    for (inf,nf) in enumerate(nfaces(p)[nf_dim(p)[end][1]])
      push!(nfacedofs, [inf])
    end
    for id in 2:length(nf_dim(p)[end])
      for (inf,nf) in enumerate(nfaces(p)[nf_dim(p)[end][id]])
        push!(nfacedofs, Int[])
      end
    end
    return nfacedofs
  end
  nodes = Gridap.Polytopes.vertices_coordinates(p)
  nfacedofs = _linear_nfacedofs(p)
  return nodes, nfacedofs
end

# Determine the high order nodes of the polytope by traversing the n-faces,
# and generate their open nodes using a reference polytope of the n-face and
# a linear reference FE on top of it
function _high_order_lagrangian_nodes_polytope(p::Polytope, order)
  vs_p = Gridap.Polytopes.vertices_coordinates(p)
  ns_float_p = [i for i in vs_p]
  ref_ps = Gridap.Polytopes.nface_ref_polytopes(p)
  # rfe_p = Gridap.RefFEs._high_order_lagrangian_reffe(p,Float64,1)
  rfe_p = LagrangianRefFE(Float64,p,1)
  nfacedofs = copy(rfe_p.nfacedofs)
  nfs = nfaces(p)
  k = length(vs_p)
  for inf_dim = 1:length(nf_dim(p)[end])-1
    nfs_dim = nf_dim(p)[end][inf_dim+1]
    nfs_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,inf_dim,0)
    for (i_nf_dim,i_nf) in enumerate(nf_dim(p)[end][inf_dim+1])
      ref_p = ref_ps[i_nf]
      # _order = Tuple(order*ones(Int,length(ref_p.extrusion)))
      _order = _extract_nonzeros(nfs[i_nf].extrusion,order)
      ns = Gridap.Polytopes._interior_nodes_int_coords(ref_p, _order)
      ns_float = Gridap.Polytopes._interior_nodes_int_to_real_coords(ns,_order)
      rfe = LagrangianRefFE(Float64,ref_p,1)
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

# Given the evaluation of linear shape functions of the reference
# polytope of a given n-face evaluated on the open high order nodes
# in the reference space, map these high order nodes in the n-face
# by combining the shape functions with the n-face vertices
function _map_high_order_lagrangian_nodes(shfs_ns, vs)
  a = shfs_ns
  T = typeof(vs[1].array)
  ndofs, npoints = size(a)
  v = Point{length(vs[1]),Float64}[]
  for j in 1:npoints
    aux = zero(T)
    for i in 1:ndofs
      aux += outer(a[i,j],vs[i]).array
    end
    cs = convert(Vector{Float64}, [aux...])
    push!(v, aux)
  end
  return v
end

# @santiagobadia : The basis for prisms, pyramids, etc, not implemented
# It would require to enrich the polynomial machinery, too compicated for a
# filter, I think it would not be efficient
function _monomial_basis(p::Polytope{D}, T, order) where D
  if (_is_hex(p))
    orders = order*ones(Int,D)
    prebasis = MonomialBasis(T,orders)
  elseif (_is_tet(p))
    filter(e,order) = sum(e) <= order
    prebasis = MonomialBasis(D, T, filter, order)
  else
    @notimplemented
  end
end

function _is_hex(p)
  for i in 2:length(extrusion(p))
    if extrusion(p)[i] != 1 return false && break end
  end
  return true
end

function _is_tet(p)
  for i in 2:length(extrusion(p))
    if extrusion(p)[i] != 2 return false && break end
  end
  return true
end

function _genereate_nface_to_dofs(nface_to_nodes,node_and_comp_to_dof)
  nnode, ncomp = size(node_and_comp_to_dof)
  nface_to_dofs = [zeros(Int,length(nodes)*ncomp) for nodes in nface_to_nodes]
  nnface = length(nface_to_nodes)
  nface_to_i = ones(Int,nnface)
  for (nface,nodes) in enumerate(nface_to_nodes)
    for comp in 1:ncomp
      for node in nodes
        dof = node_and_comp_to_dof[node,comp]
        i = nface_to_i[nface]
        nface_to_dofs[nface][i] = dof
        nface_to_i[nface] += 1
      end
    end
  end
  nface_to_dofs
end

function _extract_nonzeros(mask,values)
  b = Int[]
  for (m,n) in zip(mask,values)
    if (m != 0)
      push!(b, n)
    end
  end
  return Tuple(b)
end

end # module RefFEs
