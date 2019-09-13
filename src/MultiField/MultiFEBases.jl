module MultiFEBases

using Gridap
using Gridap.Helpers

import Gridap: gradient
import Gridap: symmetric_gradient
import Base: div
import Gridap: trace
import Gridap: curl
import Gridap: inner
import Gridap: varinner
import Base: +, -, *
import Base: length, getindex
import Gridap.FESpaces: FEBasis
import Gridap: CellBasis
import Gridap: restrict
import Gridap: Triangulation

struct FEBasisWithFieldId{B<:FEBasis}
  febasis::B
  fieldid::Int
end

Triangulation(a::FEBasisWithFieldId) = Triangulation(a.febasis)

function CellBasis(
  trian::Triangulation{D,Z},
  fun::Function,
  b::FEBasisWithFieldId,
  u...) where {D,Z}

  febasis = CellBasis(trian,fun,b.febasis,u...)
  FEBasisWithFieldId(febasis,b.fieldid)
end

for op in (:+,:-,:(gradient),:(symmetric_gradient),:(div),:(trace),:(curl))
  @eval begin
    function ($op)(a::FEBasisWithFieldId)
      FEBasisWithFieldId($op(a.febasis),a.fieldid)
    end
  end
end

for op in (:+, :-, :*)
  @eval begin

    function ($op)(a::FEBasisWithFieldId,b::CellField)
      FEBasisWithFieldId($op(a.febasis,b),a.fieldid)
    end

    function ($op)(a::CellField,b::FEBasisWithFieldId)
      FEBasisWithFieldId($op(a,b.febasis),b.fieldid)
    end

    function ($op)(a::FEBasisWithFieldId,b::Number)
      FEBasisWithFieldId($op(a.febasis,b),a.fieldid)
    end

    function ($op)(a::Number,b::FEBasisWithFieldId)
      FEBasisWithFieldId($op(a,b.febasis),b.fieldid)
    end

  end
end

function inner(a::FEBasisWithFieldId,b::CellField)
  block = inner(a.febasis,b)
  blocks = [block,]
  fieldids = [(a.fieldid,),]
  MultiCellMap(blocks,fieldids)
end

function varinner(a::FEBasisWithFieldId,b::CellField)
  inner(a,b)
end

function inner(a::FEBasisWithFieldId,f::Function)
  b = CellField(a.febasis.trian,f)
  inner(a,b)
end

function inner(a::FEBasisWithFieldId,b::FEBasisWithFieldId)
  block = inner(a.febasis,b.febasis)
  blocks = [block,]
  fieldids = [(a.fieldid,b.fieldid),]
  MultiCellMap(blocks,fieldids)
end

function varinner(a::FEBasisWithFieldId,b::FEBasisWithFieldId)
  inner(a,b)
end

function restrict(feb::FEBasisWithFieldId,trian::SkeletonTriangulation)
  sp = restrict(feb.febasis,trian)
  _new_sp(sp,feb.fieldid)
end

function _new_sp(sp::SkeletonPair{Z,T,N},fieldid) where {Z,T,N}
  b1 = FEBasisWithFieldId(sp.cellfield1,fieldid)
  b2 = FEBasisWithFieldId(sp.cellfield2,fieldid)
  SkeletonPair{Z,T,N}(b1,b2)
end

struct MultiFEBasis
  blocks::Vector{<:FEBasisWithFieldId}
end

function FEBasis(mfes::MultiFESpace)
  blocks = [ FEBasisWithFieldId(FEBasis(v),i) for (i,v) in enumerate(mfes) ]
  MultiFEBasis(blocks)
end

length(mfb::MultiFEBasis) = length(mfb.blocks)

getindex(mfb::MultiFEBasis,i::Integer) = mfb.blocks[i]

function restrict(mfeb::MultiFEBasis,trian::BoundaryTriangulation)
  blocks = [
    FEBasisWithFieldId(restrict(feb.febasis,trian),feb.fieldid)
    for feb in mfeb.blocks ]
  MultiFEBasis(blocks)
end

function restrict(mfeb::MultiFEBasis,trian::SkeletonTriangulation)
  [ restrict(feb,trian) for feb in mfeb.blocks ]
end

end # module MultiFEBases
