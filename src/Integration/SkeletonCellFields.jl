module SkeletonCellFields

using Gridap

export jump
import Gridap: restrict
import Gridap: gradient
import Gridap: varinner
import Gridap: integrate

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

struct SkeletonPair{Z,T,N}
  cellfield1::CellFieldLike{Z,T,N}
  cellfield2::CellFieldLike{Z,T,N}
end

function gradient(sp::SkeletonPair)
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  SkeletonPair(gradient(cf1),gradient(cf2))
end

function jump(sp::SkeletonPair{Z,T,1}) where {Z,T}
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  cf1 - cf2
end

function jump(sp::SkeletonPair{Z,T,2}) where {Z,T}
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  SkeletonCellBasis(cf1, -cf2)
end

struct SkeletonCellBasis{Z,T}
  cellbasis1::CellBasis{Z,T}
  cellbasis2::CellBasis{Z,T}
end

function varinner(a::SkeletonCellBasis{Z},b::CellField{Z}) where Z
  cb1 = a.cellbasis1
  cb2 = a.cellbasis2
  cm1 = varinner(cb1,b)
  cm2 = varinner(cb2,b)
  VarinnerSkeletonVector(cm1,cm2)
end

function varinner(a::SkeletonCellBasis{Z},b::SkeletonCellBasis{Z}) where Z
  a1 = a.cellbasis1
  a2 = a.cellbasis2
  b1 = b.cellbasis1
  b2 = b.cellbasis2
  cm11 = varinner(a1,b1)
  cm12 = varinner(a1,b2)
  cm21 = varinner(a2,b1)
  cm22 = varinner(a2,b2)
  VarinnerSkeletonMatrix(cm11,cm12,cm21,cm22)
end

struct VarinnerSkeletonVector{D,T}
  cellmap1::CellMap{Point{D},1,T,2}
  cellmap2::CellMap{Point{D},1,T,2}
end

struct VarinnerSkeletonMatrix{D,T}
  cellmap11::CellMap{Point{D},1,T,3}
  cellmap12::CellMap{Point{D},1,T,3}
  cellmap21::CellMap{Point{D},1,T,3}
  cellmap22::CellMap{Point{D},1,T,3}
end

function integrate(
  v::VarinnerSkeletonVector{Z},
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where Z

  v1 = integrate(v.cellmap1,trian,quad)
  v2 = integrate(v.cellmap2,trian,quad)
  VarinnerSkeletonVectorValue(v1,v2)
end



end # module
