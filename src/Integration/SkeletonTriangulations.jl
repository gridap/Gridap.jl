module SkeletonTriangulations

using Gridap

export SkeletonTriangulation

import Gridap: CellPoints
import Gridap: CellBasis
import Gridap: CellGeomap
import Gridap: CellRefFEs

struct SkeletonTriangulation{Z,D} <: Triangulation{Z,D}
  trian::Triangulation{Z,D}
  descriptor1::BoundaryDescriptor
  descriptor2::BoundaryDescriptor
  
  function SkeletonTriangulation(
    trian::Triangulation{Z,D},
    descriptor1::BoundaryDescriptor,
    descriptor2::BoundaryDescriptor) where {Z,D}
    @assert D == Z + 1
    new{Z,D}(trian,descriptor1,descriptor2)
  end

end

for op in (:CellPoints,:CellBasis,:CellGeomap,:CellRefFEs)
  @eval begin
    $op(t::SkeletonTriangulation) = $op(t.trian)
  end
end

end # module
