module BoundaryCellQuadratures

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

export BoundaryCellQuadrature
export coordinates_in_ref_cells
import Gridap: coordinates
import Gridap: weights

struct BoundaryCellQuadrature{D,Q} <: IndexCellValue{Q,1}
  quad::CellQuadrature{D}
  cell_to_polytope::IndexCellValue{<:Polytope}
  facet_to_cell::IndexCellNumber
  facet_to_lfacet::IndexCellNumber
  #TODO
  #facet_to_perm::IndexCellNumber
end

function BoundaryCellQuadrature(
  quad::IndexCellValue{<:Quadrature{D}},
  cell_to_polytope::IndexCellValue{<:Polytope},
  facet_to_cell::IndexCellNumber{<:Integer},
  facet_to_lfacet::IndexCellNumber{<:Integer}) where D
  Q = eltype(quad)
  BoundaryCellQuadrature{D,Q}(
    quad,cell_to_polytope,facet_to_cell,facet_to_lfacet)
end

coordinates(q::BoundaryCellQuadrature) = coordinates(q.quad)

weights(q::BoundaryCellQuadrature) = weights(q.quad)

function coordinates_in_ref_cells(quad::BoundaryCellQuadrature{D}) where D
 facet_qpoint_to_qcoord  = coordinates(quad)
 _coordinates_in_ref_cells(
   D,facet_qpoint_to_qcoord,quad.facet_to_lfacet,quad.cell_to_polytope)
end

function _coordinates_in_ref_cells(
  D,
  facet_qpoint_to_qcoord,facet_to_lfacet,cell_to_polytope)
  @notimplemented
end

function _coordinates_in_ref_cells(
  D,
  facet_qpoint_to_qcoord::ConstantCellValue,
  facet_to_lfacet,
  cell_to_polytope::ConstantCellValue)

  polytope = cell_to_polytope.value
  grid = Grid(polytope,D)
  trian = Triangulation(grid)
  phi = CellGeomap(trian)
  nlfacets = length(phi)

  _qpoint_to_qcoord = facet_qpoint_to_qcoord.value
  qpoint_to_qcoord = ConstantCellValue(_qpoint_to_qcoord,nlfacets)
  _lfacet_qpoint_to_rcoord = evaluate(phi,qpoint_to_qcoord)

  lfacet_qpoint_to_rcoord = CellValueFromArray(collect(_lfacet_qpoint_to_rcoord))
  facet_qpoint_to_rcoord = reindex(lfacet_qpoint_to_rcoord,facet_to_lfacet)

  facet_qpoint_to_rcoord
end

end # module
