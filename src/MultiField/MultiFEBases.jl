module MultiFEBases

using Gridap

import Gridap: gradient
import Gridap: inner
import Base: +, -, *
import Base: length, getindex
import Gridap.FESpaces: FEBasis
import Gridap: CellBasis

struct FEBasisWithFieldId{B<:CellBasis}
  cellbasis::B
  fieldid::Int
end

function CellBasis(
  trian::Triangulation{D,Z},
  fun::Function,
  b::FEBasisWithFieldId,
  u::Vararg{<:CellField{Z}}) where {D,Z}
  basis = CellBasis(trian,fun,b.cellbasis,u...)
  FEBasisWithFieldId(basis,b.fieldid)
end

for op in (:+, :-, :(gradient))
  @eval begin
    function ($op)(a::FEBasisWithFieldId)
      FEBasisWithFieldId($op(a.cellbasis),a.fieldid)
    end
  end
end

for op in (:+, :-, :*)
  @eval begin
    function ($op)(a::FEBasisWithFieldId,b::CellMap)
      FEBasisWithFieldId($op(a.cellbasis,b),a.fieldid)
    end
    function ($op)(a::CellMap,b::FEBasisWithFieldId)
      FEBasisWithFieldId($op(a,b.cellbasis),b.fieldid)
    end
  end
end

function inner(a::FEBasisWithFieldId,b::CellField)
  block = varinner(a.cellbasis,b)
  blocks = [block,]
  fieldids = [(a.fieldid,),]
  MultiCellMap(blocks,fieldids)
end

function inner(a::FEBasisWithFieldId,b::FEBasisWithFieldId)
  block = varinner(a.cellbasis,b.cellbasis)
  blocks = [block,]
  fieldids = [(a.fieldid,b.fieldid),]
  MultiCellMap(blocks,fieldids)
end

struct MultiFEBasis
  blocks::Vector{<:FEBasisWithFieldId}
end

function FEBasis(mfes::MultiFESpace)
  blocks = [ FEBasisWithFieldId(CellBasis(v),i) for (i,v) in enumerate(mfes) ]
  MultiFEBasis(blocks)
end

length(mfb::MultiFEBasis) = length(mfb.blocks)

getindex(mfb::MultiFEBasis,i::Integer) = mfb.blocks[i]

end # module MultiFEBases
