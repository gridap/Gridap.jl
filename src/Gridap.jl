"""
Gridap, grid-based approximation of PDEs in the Julia programming language

This module provides rich set of tools for the numerical solution of PDE, mainly based
on finite element methods.

The module is structured in the following sub-modules:

- [`Gridap.Helpers`](@ref)
- [`Gridap.Inference`](@ref)
- [`Gridap.Io`](@ref)
- [`Gridap.Algebra`](@ref)
- [`Gridap.TensorValues`](@ref)
- [`Gridap.Arrays`](@ref)
- [`Gridap.Fields`](@ref)
- [`Gridap.Polynomials`](@ref)
- [`Gridap.Integration`](@ref)
- [`Gridap.ReferenceFEs`](@ref)
- [`Gridap.Geometry`](@ref)
- [`Gridap.FESpaces`](@ref)
- [`Gridap.MultiField`](@ref)
- [`Gridap.Visualization`](@ref)

The exported names are:
$(EXPORTS)
"""
module Gridap

using DocStringExtensions

include("Helpers/Helpers.jl")

include("Inference/Inference.jl")

include("Io/Io.jl")

include("Algebra/Algebra.jl")

include("Arrays/Arrays.jl")

include("TensorValues/TensorValues.jl")

include("Fields/Fields.jl")

include("Polynomials/Polynomials.jl")

include("Integration/Integration.jl")

include("ReferenceFEs/ReferenceFEs.jl")

#include("Geometry/Geometry.jl")
#
#include("FESpaces/FESpaces.jl")
#
#include("MultiField/MultiField.jl")
#
#include("Visualization/Visualization.jl")
#
#include("Exports.jl")

end # module
