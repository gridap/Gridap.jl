var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Numa.jl",
    "title": "Numa.jl",
    "category": "page",
    "text": ""
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
    "location": "#Numa.FESpace",
    "page": "Numa.jl",
    "title": "Numa.FESpace",
    "category": "type",
    "text": "FE Space structure, where only one RefFE is possible in the whole mesh (to be improved in the future)\n\n\n\n\n\n"
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
