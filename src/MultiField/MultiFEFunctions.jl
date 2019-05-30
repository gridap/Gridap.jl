module MultiFEFunctions

using Gridap
using Gridap.Helpers
using Gridap.FESpaces
using Gridap.MultiAssemblers

export MultiFEFunction

import Base: length
import Base: getindex

const MultiFESpace = Vector{<:FESpace}

struct MultiFEFunction
  fields::Vector{<:FEFunction}
  free_dofs_all_fields::AbstractVector
end

function MultiFEFunction(
  free_dofs_all_fields::AbstractVector
  fespaces::Vector{<:FESpaceWithDirichletData},
  assem::MultiAssembler)
  fields = [
    FEFunction(U,restrict_to_field(assem,free_dofs_all_fields,i))
    for (i,U) in enumerate(fespaces) ]
end

length(self::MultiFEFunction) = length(self.fields)

getindex(self::MultiFEFunction,fieldid::Integer) = self.fields[fieldid]

end # module MultiFEFunctions
