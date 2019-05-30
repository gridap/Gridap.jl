module MultiFESpaces

using Gridap
using Gridap.FESpaces

export MultiFESpace

import Base: length
import Base: getindex
import Base: iterate
import Gridap.FESpaces: num_free_dofs

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

end # module MultiFESpaces
