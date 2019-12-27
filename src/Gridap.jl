"""
Gridap, grid-based approximation of PDEs in the Julia programming language

This module provides rich set of tools for the numerical solution of PDE, mainly based
on finite element methods.

The module is structured in the following sub-modules:

- [`Gridap.Helpers`](@ref)
- [`Gridap.Inference`](@ref)
- [`Gridap.Io`](@ref)
- [`Gridap.TensorValues`](@ref)
- [`Gridap.Arrays`](@ref)
- [`Gridap.Fields`](@ref)
- [`Gridap.Polynomials`](@ref)
- [`Gridap.ReferenceFEs`](@ref)
- [`Gridap.Geometry`](@ref)
- [`Gridap.Visualization`](@ref)

The exported names are:
$(EXPORTS)
"""
module Gridap

include("new/Gridap.jl")

end # module

#__precompile__()
#
#module Gridap
#
#using Reexport
#
#include("Utils/files.jl")
#
#include("CellValues/files.jl")
#
#include("Fields/files.jl")
#
#include("RefFEs/files.jl")
#
#include("Integration/files.jl")
#
#include("Geometry/files.jl")
#
#include("Algebra/files.jl")
#
#include("FESpaces/files.jl")
#
#include("MultiField/files.jl")
#
#include("Visualization/files.jl")
#
#end # module
