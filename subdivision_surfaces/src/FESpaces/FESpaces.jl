module FESpaces

using Gridap.Helpers: @check, findfirstvalue
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces: CellFE

using GridapSubdivisionSurfaces.ReferenceFEs: Loop

import Gridap.FESpaces: FESpace
import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!

include("LoopFESpaces.jl")

end
