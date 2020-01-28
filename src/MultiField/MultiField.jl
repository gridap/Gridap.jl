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
import Gridap.FESpaces: FECellBasisStyle
import Gridap.FESpaces: FEFunctionStyle
import Gridap.FESpaces: num_free_dofs
import Gridap.FESpaces: get_cell_basis
import Gridap.FESpaces: FEFunction
import Gridap.FESpaces: zero_free_values
import Gridap.FESpaces: apply_constraints_matrix_cols
import Gridap.FESpaces: apply_constraints_matrix_rows
import Gridap.FESpaces: apply_constraints_vector

export num_fields
export compute_field_offsets
export restrict_to_field

include("MultiCellArrays.jl")

include("MultiFieldCellBases.jl")

include("MultiFieldFESpaces.jl")

include("MultiFieldFEFunctions.jl")

end # module
