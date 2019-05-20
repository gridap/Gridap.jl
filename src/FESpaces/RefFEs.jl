module RefFEs

using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Maps
using Gridap.Polytopes
using Gridap.Polynomials
using Gridap.Maps: Basis
using Gridap.Maps: return_size

export DOFBasis
export RefFE
export LagrangianRefFE
export shfbasis

import Gridap: evaluate, evaluate!

# Abstract types and interfaces

"""
Abstract DOF cell basis
"""
abstract type DOFBasis{D,T} end

"""
Evaluate the DOFs for a given polynomial basis
"""
function evaluate(this::DOFBasis{D,T},
  fields::Map{Point{D},N,T,M})::Array{Float64,M} where {D,T,N,M}
  size_arr = return_size(fields, (numlocaldofs(this),))
  aux = zeros(Float64, size_arr...)
  evaluate!(this,fields,aux)
  return aux
end

function evaluate!(this::DOFBasis{D,T},
  fields::Map{Point{D},N,T,N}, auxv::AbstractArray)::Array{Float64,N} where {D,T,N}
  @abstractmethod
end
# @santiagobadia : To replace the following ones

function evaluate!(this::DOFBasis{D,T},
  prebasis::Basis{D,T}, auxv::AbstractMatrix)::Array{Float64,2} where {D,T}
  @abstractmethod
end

function evaluate!(this::DOFBasis{D,T},
  prebasis::Field{D,T}, auxv::AbstractVector)::Vector{Float64} where {D,T}
  @abstractmethod end

"""
Lagrangian DOF basis, which consists on evaluating the polynomial basis
(prebasis) in a set of points (nodes)
"""
struct LagrangianDOFBasis{D,T} <: DOFBasis{D,T}
  nodes::Array{Point{D}}
end

"""
Evaluate the Lagrangian DOFs basis (i.e., nodal values) for a given polynomial
basis
"""
function evaluate!(this::LagrangianDOFBasis{D,T},
  prebasis::Basis{D,T}, b::AbstractMatrix{Float64}) where {D,T}
  vals = Polynomials.evaluate(prebasis,this.nodes)
  l = length(prebasis); lt = length(T)
  # E = eltype(T)
  # b = Array{E,2}(undef,l, l)
  nnd = length(this.nodes)
  @assert nnd*length(T) == length(prebasis)
  function computeb!(a,b,lt,nnd)
    for k in 1:lt
      off = nnd*(k-1)
      for j in 1:size(a,2)
        for i in 1:size(a,1)
          b[i,j+off] = a[i,j][k]
        end
      end
    end
  end
  computeb!(vals,b,lt,nnd)
  return b
end

"""
Evaluate the Lagrangian DOFs basis (i.e., nodal values) for a given field in
the reference space
"""
function numlocaldofs(this::LagrangianDOFBasis{D,T}) where {D,T}
  lt = length(T)
  # E = eltype(T)
  nnd = length(this.nodes)
  # return b = Vector{E}(undef,lt*nnd)
  return lt*nnd
end

# @santiagobadia : Be careful, a physical field must be composed with geomap
# before being used here. Is this what we want?
function evaluate!(this::LagrangianDOFBasis{D,T},
  field::Field{D,T}, b::AbstractVector{Float64}) where {D,T}
  vals = Maps.evaluate(field,this.nodes)
  # I would like to use evaluate everywhere, putting evaluate in Gridap and
  # importing it in all submodules
  # This way we could use the same evaluate for bases and fields...
  # @santiagobadia : TO BE DONE
  lt = length(T)
  # E = eltype(T)
  nnd = length(this.nodes)
  # b = Vector{E}(undef,lt*nnd)
  function computeb!(a,b,lt,nnd)
    for k in 1:lt
      off = nnd*(k-1)
      for j in 1:nnd
        b[j+off] = a[j][k]
      end
    end
  end
  computeb!(vals,b,lt,nnd)
  return b
end

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
  shfbasis::PolynomialBasisWithChangeOfBasis{D,T}
  nfacedofs::Vector{Vector{Int}}
end

function LagrangianRefFE{D,T}(polytope::Polytope{D},
  orders::Array{Int64,1}) where {D,T}
  nodes=NodesArray(polytope,orders)
  dofsb = LagrangianDOFBasis{D,T}(nodes.coordinates)
  prebasis = TensorProductMonomialBasis{D,T}(orders)
  aux = zeros(Float64,numlocaldofs(dofsb),numlocaldofs(dofsb))
  @assert numlocaldofs(dofsb) == length(prebasis)
  changeofbasis=inv(evaluate!(dofsb,prebasis,aux))
  basis = PolynomialBasisWithChangeOfBasis{D,T}(prebasis, changeofbasis)
  nfacedofs=nodes.nfacenodes
  LagrangianRefFE{D,T}(polytope, dofsb, basis, nfacedofs)
end

dofs(this::LagrangianRefFE{D,T} where {D,T}) = this.dofbasis

polytope(this::LagrangianRefFE{D,T} where {D,T}) = this.polytope

shfbasis(this::LagrangianRefFE{D,T} where {D,T}) = this.shfbasis

nfacedofs(this::LagrangianRefFE{D,T} where {D,T}) = this.nfacedofs

end # module RefFEs
