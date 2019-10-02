module MultiFESpaces

using Gridap
using Gridap.Helpers
using Gridap.CachedSubVectors

export MultiFESpace

export MultiFieldStyle
export ConsequtiveMultiFieldStyle
export StridedMultiFieldStyle
export restrict_to_field

import Base: length
import Base: getindex
import Base: iterate
import Gridap.FESpaces: num_free_dofs
import Gridap.FESpaces: apply_constraints
import Gridap.FESpaces: apply_constraints_rows
import Gridap.FESpaces: apply_constraints_cols
import Gridap.FESpaces: celldofids

abstract type MultiFieldStyle end

struct ConsequtiveMultiFieldStyle <: MultiFieldStyle end

struct StridedMultiFieldStyle <: MultiFieldStyle end

struct MultiFESpace{E,S<:MultiFieldStyle}
  fespaces::Vector{<:FESpaceWithDirichletData}
end

function MultiFESpace(fespaces::Vector{<:FESpaceWithDirichletData})
  MultiFESpace(fespaces,ConsequtiveMultiFieldStyle())
end

function MultiFESpace(
  fespaces::Vector{<:FESpaceWithDirichletData}, mfs::S) where S <:MultiFieldStyle
  @assert length(fespaces) > 0
  E = eltype(value_type(fespaces[1]))
  @assert all( ( eltype(value_type(U)) == E for U in fespaces ) )
  MultiFESpace{E,S}(fespaces)
end

MultiFieldStyle(::Type{MultiFESpace{E,S}}) where {E,S} = S

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

function apply_constraints(self::MultiFESpace, mcv::MultiCellVector{T}, cellids::CellNumber) where T
  blocks = mcv.blocks
  i_to_fieldid = mcv.fieldids
  newblocks = CellVector{T}[]
  for (i,block) in enumerate(blocks)
    ifield, = i_to_fieldid[i]
    Ui = self[ifield]
    newblock = apply_constraints(Ui,block,cellids)
    push!(newblocks,newblock)
  end
  MultiCellVector(newblocks,i_to_fieldid)
end

function apply_constraints_rows(self::MultiFESpace, mcm::MultiCellMatrix{T}, cellids::CellNumber) where T
  blocks = mcm.blocks
  i_to_fieldid = mcm.fieldids
  newblocks = CellMatrix{T}[]
  for (i,block) in enumerate(blocks)
    ifield, = i_to_fieldid[i]
    Ui = self[ifield]
    newblock = apply_constraints_rows(Ui,block,cellids)
    push!(newblocks,newblock)
  end
  MultiCellMatrix(newblocks,i_to_fieldid)
end

function apply_constraints_cols(self::MultiFESpace, mcm::MultiCellMatrix{T}, cellids::CellNumber) where T
  blocks = mcm.blocks
  i_to_fieldid = mcm.fieldids
  newblocks = CellMatrix{T}[]
  for (i,block) in enumerate(blocks)
    _, jfield = i_to_fieldid[i]
    Ui = self[jfield]
    newblock = apply_constraints_cols(Ui,block,cellids)
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

function restrict_to_field(
  this::MultiFESpace{E}, v::AbstractVector{E},field::Integer) where E
  @notimplemented
end

function restrict_to_field(
  this::MultiFESpace{E,ConsequtiveMultiFieldStyle},
  v::AbstractVector{E},field::Integer) where E
  U = this.fespaces
  _restrict_to_field(U,v,field)
end

function _restrict_to_field(U,v,field)
  offsets = _compute_offsets(U)
  pini = offsets[field] + 1
  pend = offsets[field] + num_free_dofs(U[field])
  CachedSubVector(v,pini,pend)
end

function _compute_offsets(U)
  n = length(U)
  offsets = zeros(Int,n)
  for i in 1:(n-1)
    Ui = U[i]
    offsets[i+1] = offsets[i] + num_free_dofs(Ui)
  end
  offsets
end

# Pretty printing

import Base: show

function show(io::IO,self::Vector{<:FESpace})
  s = "Vector of FESpace objects with $(length(self)) fields"
  print(io,s)
end

function show(io::IO,::MIME"text/plain",self::Vector{<:FESpace})
  show(io,self)
  print(io,":")
  dofs = 0
  for U in self
    dofs += num_free_dofs(U)
  end
  print(io,"\nnfree (all fields): $dofs")
  for (i,U) in enumerate(self)
    print(io,"\nfield $i:\n")
    show(io,"text/plain",U)
  end
end

function show(io::IO,self::MultiFESpace)
  s = "MultiFESpace object with $(length(self)) fields"
  print(io,s)
end

function show(io::IO,::MIME"text/plain",self::MultiFESpace)
  show(io,self)
  print(io,":")
  for (i,U) in enumerate(self.fespaces)
    print(io,"\nfield $i:\n")
    show(io,"text/plain",U)
  end
end

end # module MultiFESpaces
