module BoundaryCellFields

using Gridap
using Gridap.Helpers
using Gridap.CellFields: NonIterableCellFieldLike

export restrict
import Base: length
import Gridap: evaluate
import Gridap: gradient

function restrict(cf::IndexCellFieldLike,desc::BoundaryDescriptor)
  BoundaryCellFieldLike(cf,desc)
end

function restrict(cf::IndexCellFieldLike{D},trian::BoundaryTriangulation{Z}) where {D,Z}
  @assert D == Z + 1
  restrict(cf,trian.descriptor)
end

struct BoundaryCellFieldLike{Z,T,N,D} <: NonIterableCellFieldLike{Z,T,N}
  cellfield::CellFieldLike{D,T,N}
  descriptor::BoundaryDescriptor
end

function BoundaryCellFieldLike(
  cellfield::IndexCellFieldLike{D,T,N}, descriptor::BoundaryDescriptor) where {D,T,N}
  Z = D - 1
  BoundaryCellFieldLike{Z,T,N,D}(cellfield,descriptor)
end

function gradient(cf::BoundaryCellFieldLike)
  cfg = gradient(cf.cellfield)
  BoundaryCellFieldLike(cfg,cf.descriptor)
end

length(bcf::BoundaryCellFieldLike) = length(bcf.descriptor.facet_to_cell)

function evaluate(bcf::BoundaryCellFieldLike{Z}, q::CellPoints{Z}) where Z

  facet_qpoint_to_qcoord = q
  facet_to_lfacet = bcf.descriptor.facet_to_lfacet
  cell_to_polytope = bcf.descriptor.cell_to_polytope
  facet_to_cell = bcf.descriptor.facet_to_cell

  r = coordinates_in_ref_cells(
    facet_qpoint_to_qcoord,
    facet_to_lfacet,
    cell_to_polytope,
    facet_to_cell)

  cf = reindex(bcf.cellfield,facet_to_cell)

  evaluate(cf,r)

end

function coordinates_in_ref_cells(
  facet_qpoint_to_qcoord,
  facet_to_lfacet,
  cell_to_polytope,
  facet_to_cell)

  _coordinates_in_ref_cells(
    facet_qpoint_to_qcoord,
    facet_to_lfacet,
    cell_to_polytope,
    facet_to_cell)

end

function  _coordinates_in_ref_cells(
  facet_qpoint_to_qcoord,
  facet_to_lfacet,
  cell_to_polytope,
  facet_to_cell)
  @notimplemented
end

function  _coordinates_in_ref_cells(
  facet_qpoint_to_qcoord::ConstantCellValue,
  facet_to_lfacet,
  cell_to_polytope::ConstantCellValue,
  facet_to_cell)

  polytope = cell_to_polytope.value
  D = length(polytope.extrusion)
  grid = Grid(polytope,D-1)
  trian = Triangulation(grid)
  phi = CellGeomap(trian)
  nlfacets = length(phi)

  _qpoint_to_qcoord = facet_qpoint_to_qcoord.value
  qpoint_to_qcoord = ConstantCellValue(_qpoint_to_qcoord,nlfacets)
  _lfacet_qpoint_to_rcoord = evaluate(phi,qpoint_to_qcoord)

  lfacet_qpoint_to_rcoord = collect(_lfacet_qpoint_to_rcoord)
  facet_qpoint_to_rcoord = CompressedCellArray(lfacet_qpoint_to_rcoord,facet_to_lfacet)

  facet_qpoint_to_rcoord
end

end # module
