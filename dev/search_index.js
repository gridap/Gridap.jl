var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Numa.jl",
    "title": "Numa.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#Numa.compose-Union{Tuple{S}, Tuple{D}, Tuple{Any,EvaluableCellArray{D,S,1}}} where S where D",
    "page": "Numa.jl",
    "title": "Numa.compose",
    "category": "method",
    "text": "Composes a lambda function f with a CellField g to provide f ∘ g. It has to be overloaded with 2 methods, one that returns the type of the result, and another one that returns the result\n\n\n\n\n\n"
},

{
    "location": "#Numa.derivative-Tuple{UnivariateMonomialBasis,Int64,Array{Float64,1}}",
    "page": "Numa.jl",
    "title": "Numa.derivative",
    "category": "method",
    "text": "Function to be eliminated in the future. Compute the numder-th derivative of a monomial at a set of 1D point\n\n\n\n\n\n"
},

{
    "location": "#Numa.evaluate!-Union{Tuple{T}, Tuple{D}, Tuple{MultivariatePolynomialBasis{D,T},Array{SArray{Tuple{D},Float64,1,D},1},Array{T,2}}} where T where D",
    "page": "Numa.jl",
    "title": "Numa.evaluate!",
    "category": "method",
    "text": "First axis of v for dofs, second for points\n\n\n\n\n\n"
},

{
    "location": "#Numa.evaluate-Tuple{Numa.EvaluableCellArrayFromBinaryOp,CellArray{StaticArrays.SArray{Tuple{D},Float64,1,D},1} where D}",
    "page": "Numa.jl",
    "title": "Numa.evaluate",
    "category": "method",
    "text": "Compute the operation op(a,b) for a EvaluableCellArrayFromBinaryOp\n\n\n\n\n\n"
},

{
    "location": "#Numa.evaluate-Tuple{Numa.EvaluableCellArrayFromComposition,CellArray{StaticArrays.SArray{Tuple{D},Float64,1,D},1} where D}",
    "page": "Numa.jl",
    "title": "Numa.evaluate",
    "category": "method",
    "text": "Evaluate a EvaluableCellArrayFromComposition\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradient-Tuple{MultivariatePolynomialBasis}",
    "page": "Numa.jl",
    "title": "Numa.gradient",
    "category": "method",
    "text": "Returns a MultivariatePolynomialBasis{TG,D} where TG is a type whose rank is one unit grater than the one of T\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradient-Tuple{Numa.EvaluableCellArray{D,T,1} where T where D}",
    "page": "Numa.jl",
    "title": "Numa.gradient",
    "category": "method",
    "text": "Returns another CellField object that represents the gradient. TG has a rank one order greater than the one of T\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradient-Tuple{Numa.EvaluableCellArray{D,T,2} where T where D}",
    "page": "Numa.jl",
    "title": "Numa.gradient",
    "category": "method",
    "text": "Returns another CellBasis object that represents the gradient TG is a value whose rank is one order grater than the one of T\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradshfsps-Tuple{LagrangianRefFE,Any}",
    "page": "Numa.jl",
    "title": "Numa.gradshfsps",
    "category": "method",
    "text": "Compute shape functions gradients (to be re-implemented)\n\n\n\n\n\n"
},

{
    "location": "#Numa.mapderivatives-Union{Tuple{T}, Tuple{D}, Tuple{EvaluableCellArray{D,T,2},EvaluableCellArray{D,SArray{Tuple{D},Float64,1,D},1}}} where T where D",
    "page": "Numa.jl",
    "title": "Numa.mapderivatives",
    "category": "method",
    "text": "Returns another CellBasis whose spatial derivatives are with respect to the coordinates of the range space of geomap\n\n\n\n\n\n"
},

{
    "location": "#Numa.shfsps-Tuple{LagrangianRefFE,Any}",
    "page": "Numa.jl",
    "title": "Numa.shfsps",
    "category": "method",
    "text": "Compute shape functions in a set of points\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellArray",
    "page": "Numa.jl",
    "title": "Numa.CellArray",
    "category": "type",
    "text": "Abstract type representing an iterable collection of Arrays{T,N}, where each array is associated to a cell.\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellBasis",
    "page": "Numa.jl",
    "title": "Numa.CellBasis",
    "category": "type",
    "text": "Abstract type that represents a cell-wise basis for a field space, where T is the type of value and D the dimension of the domain\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellBasisValues",
    "page": "Numa.jl",
    "title": "Numa.CellBasisValues",
    "category": "type",
    "text": "Abstract type that represents a function basis with value of type T evaluated at a collection of points in each cell\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellField",
    "page": "Numa.jl",
    "title": "Numa.CellField",
    "category": "type",
    "text": "Abstract type that represents a cell-wise field, where T stands for the type that represents the field at a point (e.g., scalar, vector, tensor) and D stands for the space dimension\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellFieldValues",
    "page": "Numa.jl",
    "title": "Numa.CellFieldValues",
    "category": "type",
    "text": "Abstract type that represents a field with value of type T evaluated at a collection of points in each cell\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellMatrices",
    "page": "Numa.jl",
    "title": "Numa.CellMatrices",
    "category": "type",
    "text": "Abstract type that represents a matrix of value T associated with a collection of points in each cell (typically the cell matrix at the quadrature points)\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellPoints",
    "page": "Numa.jl",
    "title": "Numa.CellPoints",
    "category": "type",
    "text": "An array of points for each cell\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellQuadrature",
    "page": "Numa.jl",
    "title": "Numa.CellQuadrature",
    "category": "type",
    "text": "Abstract type representing a collection of quadratures, one for each cell\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellScalars",
    "page": "Numa.jl",
    "title": "Numa.CellScalars",
    "category": "type",
    "text": "Abstract type that represents a scalar of value T associated with a collection of points in each cell\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellValues",
    "page": "Numa.jl",
    "title": "Numa.CellValues",
    "category": "type",
    "text": "An array of values for each cell\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellVectors",
    "page": "Numa.jl",
    "title": "Numa.CellVectors",
    "category": "type",
    "text": "Abstract type that represents a vector of value T associated with a collection of points in each cell (typically the cell rhs vector at the quadrature points)\n\n\n\n\n\n"
},

{
    "location": "#Numa.ConstantCellArray",
    "page": "Numa.jl",
    "title": "Numa.ConstantCellArray",
    "category": "type",
    "text": "Concrete implementation of CellArray, where the same array is associated to all cells. Typically, this is useful for discretizations with a single cell type.\n\n\n\n\n\n"
},

{
    "location": "#Numa.ConstantCellQuadrature",
    "page": "Numa.jl",
    "title": "Numa.ConstantCellQuadrature",
    "category": "type",
    "text": "A concrete implementation of CellQuadrature for the particular case that all cells have the same quadrature\n\n\n\n\n\n"
},

{
    "location": "#Numa.FESpace",
    "page": "Numa.jl",
    "title": "Numa.FESpace",
    "category": "type",
    "text": "FE Space structure, where only one RefFE is possible in the whole mesh (to be improved in the future)\n\n\n\n\n\n"
},

{
    "location": "#Numa.IndexableCellArray",
    "page": "Numa.jl",
    "title": "Numa.IndexableCellArray",
    "category": "type",
    "text": "Abstract type representing an indexable CellArray. By implementing a concrete IndexableCellArray, one automatically gets a type that is also iterable\n\n\n\n\n\n"
},

{
    "location": "#Numa.IntegrationDomain",
    "page": "Numa.jl",
    "title": "Numa.IntegrationDomain",
    "category": "type",
    "text": "This is the very minimum needed to describe the domain for numerical integration\n\n\n\n\n\n"
},

{
    "location": "#Numa.IntegrationMesh",
    "page": "Numa.jl",
    "title": "Numa.IntegrationMesh",
    "category": "type",
    "text": "Minimal interface for a mesh used for numerical integration\n\n\n\n\n\n"
},

{
    "location": "#Numa.LagrangianRefFE",
    "page": "Numa.jl",
    "title": "Numa.LagrangianRefFE",
    "category": "type",
    "text": "Reference Finite Element a la Ciarlet, i.e., it relies on a local function (polynomial) space, an array of nodes (DOFs), and a polytope (cell topology). The rank of the approximating field can be arbitrary. The current implementation relies on the prebasis (e.g., monomial basis of polynomials) and a change-of-basis (using the node array) to generate the canonical basis, i.e., the shape functions.\n\n\n\n\n\n"
},

{
    "location": "#Numa.MPoint",
    "page": "Numa.jl",
    "title": "Numa.MPoint",
    "category": "type",
    "text": "The mutable version of Point{D}\n\n\n\n\n\n"
},

{
    "location": "#Numa.MTensorValue",
    "page": "Numa.jl",
    "title": "Numa.MTensorValue",
    "category": "type",
    "text": "Mutable version of TensorValue{D,DD}\n\n\n\n\n\n"
},

{
    "location": "#Numa.MVectorValue",
    "page": "Numa.jl",
    "title": "Numa.MVectorValue",
    "category": "type",
    "text": "Mutable version of VectorValue{D}\n\n\n\n\n\n"
},

{
    "location": "#Numa.MultivariatePolynomialBasis",
    "page": "Numa.jl",
    "title": "Numa.MultivariatePolynomialBasis",
    "category": "type",
    "text": "Abstract type representing a multivariate polynomial basis with value of type T in a coordinate space of D dimensions\n\n\n\n\n\n"
},

{
    "location": "#Numa.NFace",
    "page": "Numa.jl",
    "title": "Numa.NFace",
    "category": "type",
    "text": "n-face of the polytope, i.e., any polytope of lower dimension (n) representing its boundary and the polytope itself (for n equal to the space dimension)\n\n\n\n\n\n"
},

{
    "location": "#Numa.NodesArray",
    "page": "Numa.jl",
    "title": "Numa.NodesArray",
    "category": "type",
    "text": "Array of nodes for a give polytope and order\n\n\n\n\n\n"
},

{
    "location": "#Numa.Point",
    "page": "Numa.jl",
    "title": "Numa.Point",
    "category": "type",
    "text": "Type representing a point of D dimensions\n\n\n\n\n\n"
},

{
    "location": "#Numa.Polytope",
    "page": "Numa.jl",
    "title": "Numa.Polytope",
    "category": "type",
    "text": "Aggregation of all n-faces that compose the polytope boundary and the polytope itself, the classification of n-faces with respect to their dimension and type\n\n\n\n\n\n"
},

{
    "location": "#Numa.Quadrature",
    "page": "Numa.jl",
    "title": "Numa.Quadrature",
    "category": "type",
    "text": "Abstract type representing a quadrature rule on a Polytope in a space of D dimensions\n\n\n\n\n\n"
},

{
    "location": "#Numa.ScalarValue",
    "page": "Numa.jl",
    "title": "Numa.ScalarValue",
    "category": "type",
    "text": "Type representing a scalar value\n\n\n\n\n\n"
},

{
    "location": "#Numa.TensorProductPolynomialBasis",
    "page": "Numa.jl",
    "title": "Numa.TensorProductPolynomialBasis",
    "category": "type",
    "text": "Multivariate polynomial basis obtained as tensor product of univariate polynomial basis per dimension\n\n\n\n\n\n"
},

{
    "location": "#Numa.TensorProductPolynomialBasis-Tuple{Array{Int64,1}}",
    "page": "Numa.jl",
    "title": "Numa.TensorProductPolynomialBasis",
    "category": "method",
    "text": "Provide a TensorProductPolynomialBasis for a vector order providing the order per dimension\n\n\n\n\n\n"
},

{
    "location": "#Numa.TensorProductQuadrature",
    "page": "Numa.jl",
    "title": "Numa.TensorProductQuadrature",
    "category": "type",
    "text": "Tensor product quadrature rule (nodes and weights) on a hyper cube [-1,1]^D\n\n\n\n\n\n"
},

{
    "location": "#Numa.TensorProductQuadratureOld",
    "page": "Numa.jl",
    "title": "Numa.TensorProductQuadratureOld",
    "category": "type",
    "text": "Tensor product quadrature rule (nodes and weights) integrating exactly 2´order´-1 polynomials\n\n\n\n\n\n"
},

{
    "location": "#Numa.TensorValue",
    "page": "Numa.jl",
    "title": "Numa.TensorValue",
    "category": "type",
    "text": "Type representing a tensor value of dimension D\n\n\n\n\n\n"
},

{
    "location": "#Numa.UnivariateMonomialBasis",
    "page": "Numa.jl",
    "title": "Numa.UnivariateMonomialBasis",
    "category": "type",
    "text": "Univariate monomial basis of a given order\n\n\n\n\n\n"
},

{
    "location": "#Numa.UnivariateMonomialBasis-Tuple{AbstractArray{Float64,1}}",
    "page": "Numa.jl",
    "title": "Numa.UnivariateMonomialBasis",
    "category": "method",
    "text": "Evaluate univariate monomial basis in a set of 1D points\n\n\n\n\n\n"
},

{
    "location": "#Numa.UnivariatePolynomialBasis",
    "page": "Numa.jl",
    "title": "Numa.UnivariatePolynomialBasis",
    "category": "type",
    "text": "Abstract basis of univariate polynomials\n\n\n\n\n\n"
},

{
    "location": "#Numa.UnivariatePolynomialBasis-Tuple{Int64}",
    "page": "Numa.jl",
    "title": "Numa.UnivariatePolynomialBasis",
    "category": "method",
    "text": "Create 1-dim univariate polynomial basis of UnivariateMonomialBasis type\n\n\n\n\n\n"
},

{
    "location": "#Numa.VectorValue",
    "page": "Numa.jl",
    "title": "Numa.VectorValue",
    "category": "type",
    "text": "Type representing a vector value of dimension D\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellArrayFromBinaryOp",
    "page": "Numa.jl",
    "title": "Numa.CellArrayFromBinaryOp",
    "category": "type",
    "text": "Concrete implementation of CellArray that represents the lazy result of applying a binary operation on two instances of CellArray the functions computevals! and computesize represents the binary operation to be performed on the two arrays. computesize(::NTuple{P,Int},NTuple{Q,Int})::NTuple{N,Int} provides the size of the result and computevals!(::Array{A,P},::Array{B,Q},::Array{T,N})::Array{T,N} computes the result. This type is essential for DRY (don\'t repeat yourself)\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellArrayFromUnaryOp",
    "page": "Numa.jl",
    "title": "Numa.CellArrayFromUnaryOp",
    "category": "type",
    "text": "Type that implements the lazy result of an unary operation on an instance of CellArray\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellBasisValuesFromSingleInterpolation",
    "page": "Numa.jl",
    "title": "Numa.CellBasisValuesFromSingleInterpolation",
    "category": "type",
    "text": "Concrete implementation for the case of the same interpolation on all cells, but arbitrary sampling points in each cell. This is typically needed for unfitted methods\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellBasisWithMappedDerivatives",
    "page": "Numa.jl",
    "title": "Numa.CellBasisWithMappedDerivatives",
    "category": "type",
    "text": "This type implements the result of mapderivatives\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellFieldValuesFromInterpolation",
    "page": "Numa.jl",
    "title": "Numa.CellFieldValuesFromInterpolation",
    "category": "type",
    "text": "General implementation Cell-wise interpolation of a field evaluated at a set of points CellFieldValues{TN} with a CellBasisValues{TB} of shape functions evaluated at the same points. Accepts any kind of CellBasisValues{TB} and CellFieldValues{TN} representing. This implementation only assumes that outer(·,·) is defined for instances of TB and TN. The result is of type T\n\n\n\n\n\n"
},

{
    "location": "#Numa.ConstantCellBasisValues",
    "page": "Numa.jl",
    "title": "Numa.ConstantCellBasisValues",
    "category": "type",
    "text": "Concrete implementation for the case of the same interpolation and the same sampling points on all cells\n\n\n\n\n\n"
},

{
    "location": "#Numa.EvaluableCellArray",
    "page": "Numa.jl",
    "title": "Numa.EvaluableCellArray",
    "category": "type",
    "text": "Abstract type that represents both a CellField and a CellBasis\n\n\n\n\n\n"
},

{
    "location": "#Numa.EvaluableCellArrayFromBinaryOp",
    "page": "Numa.jl",
    "title": "Numa.EvaluableCellArrayFromBinaryOp",
    "category": "type",
    "text": "Implements the results of a binary operation op between two instances a and b of EvaluableCellArray\n\n\n\n\n\n"
},

{
    "location": "#Numa.EvaluableCellArrayFromComposition",
    "page": "Numa.jl",
    "title": "Numa.EvaluableCellArrayFromComposition",
    "category": "type",
    "text": "Implements the composition a ∘ b of two instances a and b of EvaluableCellArray\n\n\n\n\n\n"
},

{
    "location": "#Numa.GradOfCellBasisMappedWithJacob",
    "page": "Numa.jl",
    "title": "Numa.GradOfCellBasisMappedWithJacob",
    "category": "type",
    "text": "Un-evaluated version of GradOfCellBasisValuesMappedWithJacob\n\n\n\n\n\n"
},

{
    "location": "#Numa.GradOfCellBasisValuesMappedWithJacob",
    "page": "Numa.jl",
    "title": "Numa.GradOfCellBasisValuesMappedWithJacob",
    "category": "type",
    "text": "Values of the gradients of a cell basis mapped with a given Jacobian J stands for the type of the values of the Jacobian, i.e., gradient(Point{D},Val(D))\n\n\n\n\n\n"
},

{
    "location": "#Numa.RefFE",
    "page": "Numa.jl",
    "title": "Numa.RefFE",
    "category": "type",
    "text": "Abstract Reference Finite Element\n\n\n\n\n\n"
},

{
    "location": "#Numa.jl-1",
    "page": "Numa.jl",
    "title": "Numa.jl",
    "category": "section",
    "text": "Documentation for the Numa libraryModules = [Numa,]\nOrder   = [:function, :type]"
},

]}
