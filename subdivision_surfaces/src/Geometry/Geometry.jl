module Geometry

using FillArrays
using Gridap.Helpers
using Gridap.Helpers: findfirstvalue
using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.ReferenceFEs

using GridapSubdivisionSurfaces.ReferenceFEs: _box_splines_222
using GridapSubdivisionSurfaces.ReferenceFEs

import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!

export loop_surface_model

include("SubdivisionSurfaceModels.jl")

end
