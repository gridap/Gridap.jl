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

using Gridap.Helpers
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry

using Gridap.ReferenceFEs: SerendipityPolytope

export writevtk
export write_vtk_file

include("Vtk.jl")

end # module
