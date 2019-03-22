var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Numa.jl",
    "title": "Numa.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#Numa.evaluate!-Union{Tuple{T}, Tuple{D}, Tuple{MultivariatePolynomialBasis{D,T},Array{SArray{Tuple{D},Float64,1,D},1},Array{T,2}}} where T where D",
    "page": "Numa.jl",
    "title": "Numa.evaluate!",
    "category": "method",
    "text": "First axis of v for dofs, second for points\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradient-Tuple{CellBasis}",
    "page": "Numa.jl",
    "title": "Numa.gradient",
    "category": "method",
    "text": "Returns another CellBasis object that represents the gradient TG is a value whose rank is one order grater than the one of T\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradient-Tuple{CellField}",
    "page": "Numa.jl",
    "title": "Numa.gradient",
    "category": "method",
    "text": "Returns another CellField object that represents the gradient. TG has a rank one order greater than the one of T\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradient-Tuple{MultivariatePolynomialBasis}",
    "page": "Numa.jl",
    "title": "Numa.gradient",
    "category": "method",
    "text": "Returns a MultivariatePolynomialBasis{TG,D} where TG is a type whose rank is one unit grater than the one of T\n\n\n\n\n\n"
},

{
    "location": "#Numa.gradshfsps-Tuple{LagrangianRefFE,Any}",
    "page": "Numa.jl",
    "title": "Numa.gradshfsps",
    "category": "method",
    "text": "Compute shape functions gradients (to be re-implemented)\n\n\n\n\n\n"
},

{
    "location": "#Numa.mapderivatives-Union{Tuple{T}, Tuple{D}, Tuple{CellBasis{D,T},CellField{D,SArray{Tuple{D},Float64,1,D}}}} where T where D",
    "page": "Numa.jl",
    "title": "Numa.mapderivatives",
    "category": "method",
    "text": "Returns another CellBasis, whose spatial derivatives are respect to the coordinates of the range space of geomap\n\n\n\n\n\n"
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
    "text": "Abstract type representing an iterable collection of Arrays{T,N}, where each array is associated with a cell.\n\n\n\n\n\n"
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
    "location": "#Numa.MVectorValue",
    "page": "Numa.jl",
    "title": "Numa.MVectorValue",
    "category": "type",
    "text": "Mutable version of VectorValue\n\n\n\n\n\n"
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
    "location": "#Numa.PolynomialBasis-Tuple{Int64}",
    "page": "Numa.jl",
    "title": "Numa.PolynomialBasis",
    "category": "method",
    "text": "Create 1-dim polynomial basis of type [Lagrangian, Monomial] for a given order and a set of nodes. The nodes can be equispaced or Chebyshev nodes of second kind and take values [Equispaced, Chebyshev]\n\n\n\n\n\n"
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
    "location": "#Numa.VectorValue",
    "page": "Numa.jl",
    "title": "Numa.VectorValue",
    "category": "type",
    "text": "Type representing a vector value of D components\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellArrayFromBinaryOp",
    "page": "Numa.jl",
    "title": "Numa.CellArrayFromBinaryOp",
    "category": "type",
    "text": "Concrete implementation of CellArray that represents the lazy result of applying a binary operation on two instances of CellArray the functions computevals! and computesize represents the binary operation two be performed on the two arrays. computesize(::NTuple{P,Int},NTuple{Q,Int})::NTuple{N,Int} provides the size of the result and computevals!(::Array{A,P},::Array{B,Q},::Array{T,N})::Array{T,N} computes the result. This type is essential for DRY (don\'t repeat yourself)\n\n\n\n\n\n"
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
