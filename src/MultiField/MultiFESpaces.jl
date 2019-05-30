module MultiFESpaces

using Gridap
using Gridap.CellValues
using Gridap.FESpaces
using Gridap.MultiCellArrays

export MultiFESpace

import Base: length
import Base: getindex
import Base: iterate
import Gridap.FESpaces: num_free_dofs
import Gridap.FESpaces: apply_constraints
import Gridap.FESpaces: apply_constraints_rows
import Gridap.FESpaces: apply_constraints_cols
import Gridap.FESpaces: celldofids

struct MultiFESpace{E}
  fespaces::Vector{<:FESpaceWithDirichletData}
end

function MultiFESpace(fespaces::Vector{<:FESpaceWithDirichletData})
  @assert length(fespaces) > 0
  E = eltype(value_type(fespaces[1]))
  @assert all( ( eltype(value_type(U)) == E for U in fespaces ) )
  MultiFESpace{E}(fespaces)
end

function num_free_dofs(self::MultiFESpace)
  n = 0
  for U in self
    n += num_free_dofs(U)
  end
  n
end

length(self::MultiFESpace) = length(self.fespaces)

getindex(self::MultiFESpace, i::Integer) = self.fespaces[i]

iterate(self::MultiFESpace) = iterate(self.fespaces)

iterate(self::MultiFESpace,state) = iterate(self.fespaces,state)

function apply_constraints(self::MultiFESpace, mcv::MultiCellVector{T}) where T
  blocks = mcv.blocks
  i_to_fieldid = mcv.fieldids
  newblocks = CellVector{T}[]
  for (i,block) in enumerate(blocks)
    ifield, = i_to_fieldid[i]
    Ui = self[ifield]
    newblock = apply_constraints(Ui,block)
    push!(newblocks,newblock)
  end
  MultiCellVector(newblocks,i_to_fieldid)
end

function apply_constraints_rows(self::MultiFESpace, mcm::MultiCellMatrix{T}) where T
  blocks = mcm.blocks
  i_to_fieldid = mcm.fieldids
  newblocks = CellMatrix{T}[]
  for (i,block) in enumerate(blocks)
    ifield, = i_to_fieldid[i]
    Ui = self[ifield]
    newblock = apply_constraints_rows(Ui,block)
    push!(newblocks,newblock)
  end
  MultiCellMatrix(newblocks,i_to_fieldid)
end

function apply_constraints_cols(self::MultiFESpace, mcm::MultiCellMatrix{T}) where T
  blocks = mcm.blocks
  i_to_fieldid = mcm.fieldids
  newblocks = CellMatrix{T}[]
  for (i,block) in enumerate(blocks)
    _, jfield = i_to_fieldid[i]
    Ui = self[jfield]
    newblock = apply_constraints_cols(Ui,block)
    push!(newblocks,newblock)
  end
  MultiCellMatrix(newblocks,i_to_fieldid)
end

function celldofids(self::MultiFESpace)
  n = length(self)
  blocks = Vector{CellVector{Int}}(undef,n)
  fielids = Vector{Tuple{Int}}(undef,n)
  for (i,Ui) in enumerate(self)
    blocks[i] = celldofids(Ui)
    fielids[i] = (i,)
  end
  MultiCellVector(blocks,fielids)
end

end # module MultiFESpaces
