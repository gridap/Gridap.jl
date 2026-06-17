module GridapSubdivisionSurfaces

using Gridap
using DocStringExtensions

include("ReferenceFEs/ReferenceFEs.jl")

include("Geometry/Geometry.jl")

include("FESpaces/FESpaces.jl")

macro publish(mod,name)
  quote
    using GridapSubdivisionSurfaces.$mod: $name; export $name
  end
end

@publish Geometry loop_surface_model

end
