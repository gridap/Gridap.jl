module RefFEs

using Numa.FieldValues
using Numa.Polytopes
using Numa.Polynomials

export RefFE
export LagrangianRefFE
export shfsps, gradshfsps

# Abstract types and interfaces

"""
Abstract Reference Finite Element
"""
abstract type RefFE end

# Concrete structs

"""
Reference Finite Element a la Ciarlet, i.e., it relies on a local function (polynomial) space, an array of nodes (DOFs), and a polytope (cell topology). The rank of the approximating field can be arbitrary. The current implementation relies on the prebasis (e.g., monomial basis of polynomials) and a change-of-basis (using the node array) to generate the canonical basis, i.e., the shape functions.
"""
struct LagrangianRefFE{D} <: RefFE
	polytope::Polytope{D}
	prebasis::TensorProductPolynomialBasis
	nodes::NodesArray{D}
	changeofbasis::Array{Float64,2}
	#dofs::Array{Function,1} (Other way to express it... as a function?)
	# fieldrank:: scalar, vector, tensor (or better field rank)
	# conformity, continuity (here or in the DOF numbering algorithm?)
	# dofsnface, dofsclosurenface
	dofsnface::Array{Array{Int64},1}
	rank::Int64
	numdof::Int64
end

# Methods

include("RefFEsMethods.jl")

end # module RefFEs
