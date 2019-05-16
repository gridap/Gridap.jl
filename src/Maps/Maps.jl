module Maps

using Gridap.Helpers
using Gridap.FieldValues
using Gridap.CachedArrays
using LinearAlgebra: inv

using Base.Cartesian: @nloops, @nexprs, @nref

export Map
export Field
export AnalyticalMap
export AnalyticalField
export Geomap
export Basis

import Gridap: evaluate, gradient, ∇
import Gridap: evaluate!, return_size
import Gridap: inner, outer, varinner, compose, lincomb, attachgeomap

import Base: +, -, *, /, ∘

include("AbstractMaps.jl")
include("Operations.jl")
include("Composition.jl")
include("AnalyticalField.jl")
include("Testers.jl")

end # module Maps
