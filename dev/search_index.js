var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Numa.jl",
    "title": "Numa.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#Numa.evaluate!-Union{Tuple{D}, Tuple{T}, Tuple{MultivariatePolynomialBasis{T,D},Array{SArray{Tuple{D},Float64,1,D},1},Array{T,2}}} where D where T",
    "page": "Numa.jl",
    "title": "Numa.evaluate!",
    "category": "method",
    "text": "First axis of v for dofs, second for points\n\n\n\n\n\n"
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
    "location": "#Numa.CellBasisValues",
    "page": "Numa.jl",
    "title": "Numa.CellBasisValues",
    "category": "type",
    "text": "Abstract type that represents a function basis with value of type T evaluated at a collection of points in each cell\n\n\n\n\n\n"
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
    "text": "Abstract type representing an indexable CellArray. By implementing a concrete IndexableCellArray, one automatically gets an type that is also iterable\n\n\n\n\n\n"
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
    "location": "#Numa.CellBasisValuesFromSingleInterpolation",
    "page": "Numa.jl",
    "title": "Numa.CellBasisValuesFromSingleInterpolation",
    "category": "type",
    "text": "Concrete implementation for the case of the same interpolation on all cells, but arbitrary sampling points in each cell. This is typically needed for unfitted methods\n\n\n\n\n\n"
},

{
    "location": "#Numa.CellPoints",
    "page": "Numa.jl",
    "title": "Numa.CellPoints",
    "category": "type",
    "text": "An array of points for each cell\n\n\n\n\n\n"
},

{
    "location": "#Numa.ConstantCellBasisValues",
    "page": "Numa.jl",
    "title": "Numa.ConstantCellBasisValues",
    "category": "type",
    "text": "Concrete implementation for the case of the same interpolation and the same sampling points on all cells\n\n\n\n\n\n"
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
