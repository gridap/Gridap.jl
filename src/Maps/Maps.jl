module Maps

using Numa.Helpers
using Numa.FieldValues
using Numa.CachedArrays
using LinearAlgebra: inv

using Base.Cartesian: @nloops, @nexprs, @nref

export Map
export Field
export AnalyticalMap
export AnalyticalField
export Geomap
export Basis
export varinner
export lincomb
export compose
export attachgeomap

import Numa: evaluate, gradient, ∇
import Numa: evaluate!, return_size
import Numa: inner, outer

import Base: +, -, *, /, ∘

include("AbstractMaps.jl")
include("Operations.jl")
include("Composition.jl")
include("AnalyticalField.jl")
include("Testers.jl")

end # module Maps
