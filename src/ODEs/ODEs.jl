"""

The exported names are
$(EXPORTS)
"""
module ODEs

using DocStringExtensions
using LinearAlgebra

include("ODETools/ODETools.jl")

include("TransientFETools/TransientFETools.jl")

# include("DiffEqsWrappers/DiffEqsWrappers.jl")

end #module

const GridapODEs = ODEs
