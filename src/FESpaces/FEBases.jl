module FEBases

using Gridap

export FEBasis

import Gridap: CellBasis
import Gridap: gradient
import Gridap: symmetric_gradient
import Base: div
import Gridap: trace
import Gridap: curl
import Gridap: inner
import Gridap: varinner
import Base: +, -, *
import Gridap: restrict
import Gridap: Triangulation

struct FEBasis{B<:CellBasis,T<:Triangulation}
  cellbasis::B
  trian::T
end

function FEBasis(fespace::FESpace)
  b = CellBasis(fespace)
  trian = Triangulation(fespace)
  FEBasis(b,trian)
end

for op in (:+,:-,:(gradient),:(symmetric_gradient),:(div),:(trace),:(curl))
  @eval begin
    function ($op)(a::FEBasis)
      FEBasis($op(a.cellbasis),a.trian)
    end
  end
end

for op in (:+, :-, :*)
  @eval begin

    function ($op)(a::FEBasis,b::CellField)
      FEBasis($op(a.cellbasis,cellnewaxis(b,dim=1)),a.trian)
    end

    function ($op)(a::CellField,b::FEBasis)
      FEBasis($op(cellnewaxis(a,dim=1),b.cellbasis),b.trian)
    end

    function ($op)(a::FEBasis,b::Number)
      FEBasis($op(a.cellbasis,b),a.trian)
    end

    function ($op)(a::Number,b::FEBasis)
      FEBasis($op(a,b.cellbasis),b.trian)
    end

  end
end

function inner(a::FEBasis,b::CellField)
  varinner(a.cellbasis,b)
end

function varinner(a::FEBasis,b::CellField)
  varinner(a.cellbasis,b)
end

function inner(a::FEBasis,f::Function)
  b = CellField(a.trian,f)
  inner(a,b)
end

function inner(a::FEBasis,b::FEBasis)
  varinner(a.cellbasis,b.cellbasis)
end

function varinner(a::FEBasis,b::FEBasis)
  varinner(a.cellbasis,b.cellbasis)
end

Triangulation(a::FEBasis) = a.trian

function CellBasis(
  trian::Triangulation{D,Z},
  fun::Function,
  b::FEBasis,
  u...) where {D,Z}
  basis = CellBasis(trian,fun,b.cellbasis,u...)
  FEBasis(basis,trian)
end

function restrict(feb::FEBasis,trian::BoundaryTriangulation)
  cb = restrict(feb.cellbasis,trian)
  FEBasis(cb,trian)
end

function restrict(feb::FEBasis,trian::SkeletonTriangulation)
  sp = restrict(feb.cellbasis,trian)
  _new_sp(sp,trian)
end

function _new_sp(sp::SkeletonPair{Z,T,N},trian) where {Z,T,N}
  b1 = FEBasis(sp.cellfield1,trian)
  b2 = FEBasis(sp.cellfield2,trian)
  SkeletonPair{Z,T,N}(b1,b2)
end

end # module
