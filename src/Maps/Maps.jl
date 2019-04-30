module Maps

using Numa.Helpers
using Numa.FieldValues
using Numa.CachedArrays

export Map
export Field
export AnalyticalField
export GeoMap
export Basis

import Numa: evaluate, gradient
import Numa: evaluate!, return_size
import Numa: inner, outer

import Base: +, -, *, /, ∘

include("AbstractMaps.jl")
include("Operations.jl")
include("Composition.jl")
include("AnalyticalField.jl")

end # module Maps
