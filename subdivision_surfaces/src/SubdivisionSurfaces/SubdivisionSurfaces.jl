module SubdivisionSurfaces

using FillArrays: Fill
using Gridap
using Gridap.Helpers: findfirstvalue, @check
using Gridap.Arrays: evaluate!, return_cache, Table
using Gridap.Fields: Field, GenericField, linear_combination
using Gridap.ReferenceFEs: TRI, LagrangianRefFE
using Gridap.CellData: Measure
using Gridap.Geometry: UnstructuredDiscreteModel, DiscreteModel, GridTopology,
  get_polytopes, get_faces, get_vertex_coordinates, get_grid_topology,
  UnstructuredGridTopology, FaceLabeling, OrientationStyle, UnstructuredGrid
using Gridap.FESpaces: FESpace

using GridapSubdivisionSurfaces.ReferenceFEs
using GridapSubdivisionSurfaces.ReferenceFEs: _box_splines_222
using GridapSubdivisionSurfaces.FESpaces: LoopPatchVerticesMap

export loop_surface_model

include("LoopSurfaces.jl")

end
