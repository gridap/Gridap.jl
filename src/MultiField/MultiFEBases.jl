module MultiFEBases

using Gridap
using Gridap.Helpers

import Gridap: gradient
import Gridap: inner
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
  u::Vararg{<:CellField{Z}}) where {D,Z}

  febasis = CellBasis(trian,fun,b.febasis,u...)
  FEBasisWithFieldId(febasis,b.fieldid)
end

for op in (:+, :-, :(gradient))
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
  end
end

function inner(a::FEBasisWithFieldId,b::CellField)
  block = inner(a.febasis,b)
  blocks = [block,]
  fieldids = [(a.fieldid,),]
  MultiCellMap(blocks,fieldids)
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

function restrict(feb::MultiFEBasis,trian::SkeletonTriangulation)
  @notimplemented
  # We still need to create a MultiSkeletonPair
end

end # module MultiFEBases
