module BoundaryTriangulations

using Gridap

export BoundaryTriangulation

import Gridap: CellPoints
import Gridap: CellBasis
import Gridap: CellGeomap
import Gridap: CellRefFEs

struct BoundaryTriangulation{Z,D} <: Triangulation{Z,D}
  trian::Triangulation{Z,D}
  descriptor::BoundaryDescriptor
  
  function BoundaryTriangulation(
    trian::Triangulation{Z,D}, descriptor::BoundaryDescriptor) where {Z,D}
    @assert D == Z + 1
    new{Z,D}(trian,descriptor)
  end

end

for op in (:CellPoints,:CellBasis,:CellGeomap,:CellRefFEs)
  @eval begin
    $op(t::BoundaryTriangulation) = $op(t.trian)
  end
end

end # module
