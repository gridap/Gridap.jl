module MultiField

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.Integration

using Gridap.FESpaces: _operate_cell_basis
using Gridap.FESpaces: _operate_cell_matrix_field
using Gridap.FESpaces: SkeletonCellBasis

import Gridap.Helpers: operate
import Gridap.Arrays: get_array
import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Geometry: get_cell_map
import Gridap.Geometry: similar_object
import Gridap.Geometry: restrict
import Gridap.Fields: integrate
import Gridap.FESpaces: TrialStyle

export num_fields

include("MultiCellArrays.jl")

include("MultiFieldCellBases.jl")

end # module
