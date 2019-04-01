module IntegrationMeshes

export IntegrationMesh
export geomap, cellcoordinates, cellbasis, ncells

using Numa.Helpers
using Numa.FieldValues
using Numa.CellFunctions

"""
Minimal interface for a mesh used for numerical integration
"""
abstract type IntegrationMesh{Z,D} end

cellcoordinates(::IntegrationMesh{Z,D} where {Z,D})::CellPoints{D} = @abstractmethod

cellbasis(::IntegrationMesh{Z,D} where {Z,D})::CellBasis{Z,Float64} = @abstractmethod

function geomap(self::IntegrationMesh)
  coords = cellcoordinates(self)
  basis = cellbasis(self)
  expand(basis,coords)
end

function ncells(self::IntegrationMesh)
  coords = cellcoordinates(self)
  length(coords)
end

end # module IntegrationMeshes
