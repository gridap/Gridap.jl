module FEBases

using Gridap

export FEBasis

import Gridap: CellBasis
import Gridap: gradient
import Gridap: symmetric_gradient
import Gridap: inner
import Base: +, -, *
import Gridap: restrict

struct FEBasis{B<:CellBasis,T<:Triangulation}
  cellbasis::B
  trian::T
end

function FEBasis(fespace::FESpace)
  b = CellBasis(fespace)
  trian = Triangulation(fespace)
  FEBasis(b,trian)
end

for op in (:+, :-, :(gradient),:(symmetric_gradient))
  @eval begin
    function ($op)(a::FEBasis)
      FEBasis($op(a.cellbasis),a.trian)
    end
  end
end

for op in (:+, :-, :*)
  @eval begin
    function ($op)(a::FEBasis,b::CellMap)
      FEBasis($op(a.cellbasis,b),a.trian)
    end
    function ($op)(a::CellMap,b::FEBasis)
      FEBasis($op(a,b.cellbasis),b.trian)
    end
  end
end

function inner(a::FEBasis,b::CellField)
  varinner(a.cellbasis,b)
end

function inner(a::FEBasis,f::Function)
  b = CellField(a.trian,f)
  inner(a,b)
end

function inner(a::FEBasis,b::FEBasis)
  varinner(a.cellbasis,b.cellbasis)
end

function CellBasis(
  trian::Triangulation{D,Z},
  fun::Function,
  b::FEBasis,
  u::Vararg{<:CellField{Z}}) where {D,Z}
  basis = CellBasis(trian,fun,b.cellbasis,u...)
  FEBasis(basis,trian)
end

function restrict(feb::FEBasis,trian::BoundaryTriangulation)
  cb = restrict(feb.cellbasis,trian)
  FEBasis(cb,trian)
end

function restrict(feb::FEBasis,trian::SkeletonTriangulation)
  restrict(feb.cellbasis,trian)
end

end # module
