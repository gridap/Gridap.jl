"""

The exported names are
$(EXPORTS)
"""
module Visualization

using DocStringExtensions
using WriteVTK
using WriteVTK.VTKCellTypes: VTK_VERTEX
using WriteVTK.VTKCellTypes: VTK_LINE
using WriteVTK.VTKCellTypes: VTK_TRIANGLE
using WriteVTK.VTKCellTypes: VTK_QUAD
using WriteVTK.VTKCellTypes: VTK_TETRA
using WriteVTK.VTKCellTypes: VTK_HEXAHEDRON
using WriteVTK.VTKCellTypes: VTK_QUADRATIC_TRIANGLE
using WriteVTK.VTKCellTypes: VTK_QUADRATIC_TETRA
using WriteVTK.VTKCellTypes: VTK_QUADRATIC_EDGE
using WriteVTK.VTKCellTypes: VTK_QUADRATIC_HEXAHEDRON
#using WriteVTK.VTKCellTypes: VTK_BIQUADRATIC_HEXAHEDRON
using WriteVTK.VTKCellTypes: VTK_QUADRATIC_QUAD
using WriteVTK.VTKCellTypes: VTK_BIQUADRATIC_QUAD
using WriteVTK.VTKCellTypes: VTK_WEDGE
using WriteVTK.VTKCellTypes: VTK_PYRAMID
using WriteVTK.VTKCellTypes: VTK_LAGRANGE_QUADRILATERAL
using WriteVTK.VTKCellTypes: VTK_LAGRANGE_TRIANGLE
using WriteVTK.VTKCellTypes: VTK_LAGRANGE_HEXAHEDRON
using WriteVTK.VTKCellTypes: VTK_LAGRANGE_TETRAHEDRON
using WriteVTK.VTKCellTypes: VTK_POLYGON, VTK_POLYHEDRON

using Gridap.Helpers
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.ReferenceFEs: _find_unique_with_indices
using Gridap.ReferenceFEs: SerendipityPolytope
#using Gridap.FESpaces
using FillArrays
using Gridap.CellData

using Gridap.Algebra: length_to_ptrs!, rewind_ptrs!

export writevtk
export createvtk
export createpvd
export savepvd
export write_vtk_file
export visualization_data
export VisualizationData

include("VisualizationData.jl")
include("Vtk.jl")

end # module
