module RefFEs

using Numa.FieldValues
using Numa.Polytopes
using Numa.Polynomials

export RefFE
export LagrangianRefFE
export shfsps, gradshfsps

# Abstract types and interfaces

"""
Abstract DOF cell basis
"""
abstract type DOFBasis{D,T} end

evaluate(self::DOFBasis{D,T}, prebasis::MultivariatePolynomialBasis{D,T})::Array{Float64,2} = @abstractmethod

# evaluate(self::DOFBasis, prebasis::Field{D,T})::Vector{Float64} = @abstractmethod
# Field to be implemented, to answer evaluate and gradient, it can be a local FE
# function or an analytical function

"""
Lagrangian DOF basis, which consists on evaluating the polynomial basis
(prebasis) in a set of points (nodes)
"""
struct LagrangianDOFBasis{D,T} <: DOFBasis{D,T}
  nodes::Array{Point{D}}
end

function evaluatedofs(this::LagrangianDOFBasis, 
	prebasis::MultivariatePolynomialBasis{D,T}) where {D,T}
  return evaluate(prebasis,this.nodes)
end

"""
Abstract Reference Finite Element
"""
abstract type RefFE{D,T} end

dofs(self::RefFE{D,T} where {D,T})::DOFBasis{D,T} = @abstractmethod

permutation(self::RefFE, nf::Int, cell_vertex_gids::AbstractVector{Int},
nface_vertex_gids::AbstractVector{Int}, nface_order::Int)::Vector{Int}
= @abstractmethod

polytope(self::RefFE{D,T} where {D,T})::Polytope{D} = @abstractmethod

shapefunctions(self::RefFE{D,T} where {D,T})::MultivariatePolynomialBasis{D,T}
= @abstractmethod

nfacetoowndofs(self::RefFE{D,T} where {D,T})::Vector{Vector{Int}}
= @abstractmethod

# Put changeofbasis inside concrete MultivariatePolynomialBasis

# Result is the permutation from the VEF to the local nface index within the cell
# cell_dofs[nface_dofs[v]] = face_dofs

# Field evaluate and gradient

# Concrete structs

"""
Reference Finite Element a la Ciarlet, i.e., it relies on a local function
(polynomial) space, an array of nodes (DOFs), and a polytope (cell topology).
The rank of the approximating field can be arbitrary. The current implementation
relies on the prebasis (e.g., monomial basis of polynomials) and a
change-of-basis (using the node array) to generate the canonical basis, i.e.,
the shape functions.
"""
struct LagrangianRefFE{D,T} <: RefFE{D,T}
	polytope::Polytope{D,T}
	basis::TensorProductPolynomialBasisWithChangeOfBasis
	dofs::LagrangianDofBasis
	nfacetoowndofs::Vector{Vector{Int}}
end

struct LagrangianDOFBasis <: DOFBasis
	nodes::NodesArray{D} # In dofbasis
	# send it to the polynomial basis...
	changeofbasis::Array{Float64,2} # In dofbasis
end

# Methods

include("RefFEsMethods.jl")

end # module RefFEs
