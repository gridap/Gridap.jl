"""
Gridap, grid-based approximation of PDEs in the Julia programming language

This module provides rich set of tools for the numerical solution of PDE, mainly based
on finite element methods.

The module is structured in the following sub-modules:

- [`Gridap.Helpers`](@ref)
- [`Gridap.Io`](@ref)
- [`Gridap.Algebra`](@ref)
- [`Gridap.Arrays`](@ref)
- [`Gridap.TensorValues`](@ref)
- [`Gridap.Fields`](@ref)
- [`Gridap.Polynomials`](@ref)
- [`Gridap.Integration`](@ref)
- [`Gridap.ReferenceFEs`](@ref)
- [`Gridap.CellData`](@ref)
- [`Gridap.Geometry`](@ref)
- [`Gridap.Visualization`](@ref)
- [`Gridap.FESpaces`](@ref)
- [`Gridap.MultiField`](@ref)

The exported names are:
$(EXPORTS)
"""
module Gridap

using DocStringExtensions

include("Helpers/Helpers.jl")

include("Io/Io.jl")

include("Algebra/Algebra.jl")

include("Arrays/Arrays.jl")

include("TensorValues/TensorValues.jl")

include("Fields/Fields.jl")

include("Polynomials/Polynomials.jl")

include("Integration/Integration.jl")

include("ReferenceFEs/ReferenceFEs.jl")

include("Geometry/Geometry.jl")

include("CellData/CellData.jl")

include("Visualization/Visualization.jl")

include("FESpacesNEW/FESpaces.jl")

include("MultiField/MultiField.jl")

include("Exports.jl")

end # module
