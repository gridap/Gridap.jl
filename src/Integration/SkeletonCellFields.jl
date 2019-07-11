module SkeletonCellFields

using Gridap

using Gridap.CellFields: NonIterableCellFieldLike

export jump
import Gridap: restrict
import Base: length
import Gridap: evaluate
import Gridap: gradient

function restrict(
  cf::IndexCellFieldLike,desc1::BoundaryDescriptor,desc2::BoundaryDescriptor)
  cf1 = restrict(cf,desc1)
  cf2 = restrict(cf,desc2)
  SkeletonPair(cf1,cf2)
end

function restrict(
  cf::IndexCellFieldLike{D},trian::SkeletonTriangulation{Z}) where {D,Z}
  @assert D == Z + 1
  restrict(cf,trian.descriptor1,trian.descriptor2)
end

struct SkeletonPair
  cellfield1::CellFieldLike
  cellfield2::CellFieldLike
end

function gradient(sp::SkeletonPair)
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  SkeletonPair(gradient(cf1),gradient(cf2))
end

function jump(sp::SkeletonPair)
  cf1 = sp.cellfield1
  cf2 = -sp.cellfield2
  SkeletonCellFieldLike(cf1,cf2)
end

struct SkeletonCellFieldLike{Z,T,N} <: NonIterableCellFieldLike{Z,T,N}
  cellfield1::CellFieldLike{Z,T,N}
  cellfield2::CellFieldLike{Z,T,N}
end

length(bcf::SkeletonCellFieldLike) = length(bcf.pair.cellfield1)

function evaluate(scf::SkeletonCellFieldLike{Z,T,1}, q::CellPoints{Z}) where {Z,T}
  cf1 = scf.cellfield1
  cf2 = scf.cellfield2
  v1 = evaluate(cf1,q)
  v2 = evaluate(cf2,q)
  v1 + v2
end

end # module
