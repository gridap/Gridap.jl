module BoundaryTriangulations

using Gridap

import Gridap: CellPoints
import Gridap: CellBasis
import Gridap: CellGeomap
import Gridap: CellRefFEs

struct BoundaryTriangulation{Z,D} <: Triangulation{Z,D}
  trian::Triangulation{Z,D}
  facet_to_cell::IndexCellNumber
  facet_to_lfacet::IndexCellNumber
  #TODO
  #facet_to_perm::IndexCellNumber
  
  function BoundaryTriangulation(
    trian::Triangulation{Z,D},
    facet_to_cell::IndexCellNumber,
    facet_to_lfacet::IndexCellNumber) where {Z,D}
    @assert D == Z + 1
    new{Z,D}(trian,facet_to_cell,facet_to_lfacet)
  end

end

for op in (:CellPoints,:CellBasis,:CellGeomap,:CellRefFEs)
  @eval begin
    $op(t::BoundaryTriangulation) = $op(t.trian)
  end
end

end # module
