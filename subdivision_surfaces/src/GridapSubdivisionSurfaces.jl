module GridapSubdivisionSurfaces

using Gridap
using DocStringExtensions

include("ReferenceFEs/ReferenceFEs.jl")

include("FESpaces/FESpaces.jl")

include("SubdivisionSurfaces/SubdivisionSurfaces.jl")


macro publish(mod,name)
  quote
    using GridapSubdivisionSurfaces.$mod: $name; export $name
  end
end

@publish SubdivisionSurfaces loop_surface_model

end
