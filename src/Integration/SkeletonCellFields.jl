module SkeletonCellFields

using Gridap

export SkeletonPair
export SkeletonCellBasis
export SkeletonCellVector
export SkeletonCellMatrix

export jump
export mean
import Gridap: restrict
import Gridap: inner
import Gridap: integrate
import Base: +, -, *
import Gridap: gradient
import Gridap: symmetric_gradient
import Base: div
import Gridap: trace
import Gridap: curl

function restrict(
  cf::IndexCellFieldLike{Z,T,N},
  desc1::BoundaryDescriptor,
  desc2::BoundaryDescriptor) where {Z,T,N}

  cf1 = restrict(cf,desc1)
  cf2 = restrict(cf,desc2)
  SkeletonPair{Z-1,T,N}(cf1,cf2)
end

function restrict(
  cf::IndexCellFieldLike{D},trian::SkeletonTriangulation{Z}) where {D,Z}
  @assert D == Z + 1
  restrict(cf,trian.descriptor1,trian.descriptor2)
end

struct SkeletonPair{Z,T,N}
  cellfield1
  cellfield2
end

for op in (:+,:-,:(gradient),:(symmetric_gradient),:(div),:(trace),:(curl))
  @eval begin
    function ($op)(a::SkeletonPair{Z,T,N}) where {Z,T,N}
      cf1 = a.cellfield1
      cf2 = a.cellfield2
      SkeletonPair{Z,T,N}($op(cf1),$op(cf2))
    end
  end
end

for op in (:+, :-, :*)
  @eval begin

    function ($op)(a::SkeletonPair{Z,T,N},b::CellField) where {Z,T,N}
      cf1 = a.cellfield1
      cf2 = a.cellfield2
      SkeletonPair{Z,T,N}($op(cf1,b),$op(cf2,b))
    end

    function ($op)(a::CellField,b::SkeletonPair{Z,T,N}) where {Z,T,N}
      cf1 = b.cellfield1
      cf2 = b.cellfield2
      SkeletonPair{Z,T,N}($op(a,cf1),$op(a,cf2))
    end

  end
end

function jump(sp::SkeletonPair{Z,T,1}) where {Z,T}
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  cf1 - cf2
end

function jump(sp::SkeletonPair{Z,T,2}) where {Z,T}
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  SkeletonCellBasis{Z,T}(cf1, -cf2)
end

function mean(sp::SkeletonPair{Z,T,1}) where {Z,T}
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  0.5*(cf1 + cf2)
end

function mean(sp::SkeletonPair{Z,T,2}) where {Z,T}
  cf1 = sp.cellfield1
  cf2 = sp.cellfield2
  SkeletonCellBasis{Z,T}(0.5*cf1, 0.5*cf2)
end

struct SkeletonCellBasis{Z,T}
  cellbasis1
  cellbasis2
end

function inner(a::SkeletonCellBasis{Z},b::CellField{Z}) where Z
  cb1 = a.cellbasis1
  cb2 = a.cellbasis2
  cm1 = varinner(cb1,b)
  cm2 = varinner(cb2,b)
  SkeletonVarinnerVector{Z}(cm1,cm2)
end

function inner(a::SkeletonCellBasis{Z},b::SkeletonCellBasis{Z}) where Z
  a1 = a.cellbasis1
  a2 = a.cellbasis2
  b1 = b.cellbasis1
  b2 = b.cellbasis2
  cm11 = varinner(a1,b1)
  cm12 = varinner(a1,b2)
  cm21 = varinner(a2,b1)
  cm22 = varinner(a2,b2)
  SkeletonVarinnerMatrix{Z}(cm11,cm12,cm21,cm22)
end

struct SkeletonVarinnerVector{D}
  cellmap1
  cellmap2
end

for op in (:+, :-)
  @eval begin

    function ($op)(a::SkeletonVarinnerVector{Z},b::SkeletonVarinnerVector{Z}) where Z
      c1 = $op(a.cellmap1,b.cellmap1)
      c2 = $op(a.cellmap2,b.cellmap2)
      SkeletonVarinnerVector{Z}(c1,c2)
    end
    
    function ($op)(a::SkeletonVarinnerVector{Z}) where Z
      c1 = $op(a.cellmap1)
      c2 = $op(a.cellmap2)
      SkeletonVarinnerVector{Z}(c1,c2)
    end

  end
end

function (*)(a::Real,b::SkeletonVarinnerVector{Z}) where Z
  c1 = a*b.cellmap1
  c2 = a*b.cellmap2
  SkeletonVarinnerVector{Z}(c1,c2)
end

struct SkeletonCellVector
  cellvector1
  cellvector2
end

struct SkeletonVarinnerMatrix{D}
  cellmap11
  cellmap12
  cellmap21
  cellmap22
end

for op in (:+, :-)
  @eval begin

    function ($op)(a::SkeletonVarinnerMatrix{Z},b::SkeletonVarinnerMatrix{Z}) where Z
      cellmap11 = $op(a.cellmap11, b.cellmap11)
      cellmap12 = $op(a.cellmap12, b.cellmap12)
      cellmap21 = $op(a.cellmap21, b.cellmap21)
      cellmap22 = $op(a.cellmap22, b.cellmap22)
      SkeletonVarinnerMatrix{Z}(
        cellmap11,
        cellmap12,
        cellmap21,
        cellmap22)
    end

  end
end

function (*)(a::Real,b::SkeletonVarinnerMatrix{Z}) where Z
      cellmap11 = a*b.cellmap11
      cellmap12 = a*b.cellmap12
      cellmap21 = a*b.cellmap21
      cellmap22 = a*b.cellmap22
      SkeletonVarinnerMatrix{Z}(
        cellmap11,
        cellmap12,
        cellmap21,
        cellmap22)
end

struct SkeletonCellMatrix
  cellmatrix11
  cellmatrix12
  cellmatrix21
  cellmatrix22
end

function integrate(
  v::SkeletonVarinnerVector{Z},
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where Z

  v1 = integrate(v.cellmap1,trian,quad)
  v2 = integrate(v.cellmap2,trian,quad)
  SkeletonCellVector(v1,v2)
end

function integrate(
  v::SkeletonVarinnerMatrix{Z},
  trian::Triangulation{Z},
  quad::CellQuadrature{Z}) where Z

  m11 = integrate(v.cellmap11,trian,quad)
  m12 = integrate(v.cellmap12,trian,quad)
  m21 = integrate(v.cellmap21,trian,quad)
  m22 = integrate(v.cellmap22,trian,quad)
  SkeletonCellMatrix(m11,m12,m21,m22)
end

end # module
