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
  orders::Array{Int64,1}) where {D,T}
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

dofs(this::LagrangianRefFE{D,T} where {D,T}) = this.dofbasis

polytope(this::LagrangianRefFE{D,T} where {D,T}) = this.polytope

shfbasis(this::LagrangianRefFE{D,T} where {D,T}) = this.shfbasis

nfacedofs(this::LagrangianRefFE{D,T} where {D,T}) = this.nfacedofs

end # module RefFEs
