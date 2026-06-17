module FESpaces

using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces: CellFE

using GridapSubdivisionSurfaces.ReferenceFEs: Loop
using GridapSubdivisionSurfaces.Geometry: LoopPatchVerticesMap

import Gridap.FESpaces: FESpace


include("LoopFESpaces.jl")

end
