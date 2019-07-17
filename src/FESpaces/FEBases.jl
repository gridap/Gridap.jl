module FEBases

using Gridap

export FEBasis

import Gridap: CellBasis
import Gridap: gradient
import Gridap: inner
import Base: +, -, *
import Gridap: restrict

struct FEBasis{B<:CellBasis}
  cellbasis::B
end

function FEBasis(fespace::FESpace)
  b = CellBasis(fespace)
  FEBasis(b)
end

for op in (:+, :-, :(gradient))
  @eval begin
    function ($op)(a::FEBasis)
      FEBasis($op(a.cellbasis))
    end
  end
end

for op in (:+, :-, :*)
  @eval begin
    function ($op)(a::FEBasis,b::CellMap)
      FEBasis($op(a.cellbasis,b))
    end
    function ($op)(a::CellMap,b::FEBasis)
      FEBasis($op(a,b.cellbasis))
    end
  end
end

function inner(a::FEBasis,b::CellField)
  varinner(a.cellbasis,b)
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
  FEBasis(basis)
end

function restrict(feb::FEBasis,trian::BoundaryTriangulation)
  cb = restrict(feb.cellbasis,trian)
  FEBasis(cb)
end

function restrict(feb::FEBasis,trian::SkeletonTriangulation)
  restrict(feb.cellbasis,trian)
end

end # module
