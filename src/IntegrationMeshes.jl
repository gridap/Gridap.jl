module IntegrationMeshes

export IntegrationMesh
export geomap, cellcoordinates, cellbasis, ncells
export integrate

using LinearAlgebra: det

using Numa.Helpers
using Numa.FieldValues
using Numa.Polynomials
using Numa.CellValues
using Numa.CellFunctions
using Numa.CellQuadratures
using Numa.Quadratures

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

function integrate(cellfun::CellFunction{Point{D},1,T,N},phi::CellGeomap{D,Z},quad::CellQuadrature{D}) where {D,Z,T,N}
  z = coordinates(quad)
  w = weights(quad)
  f = evaluate(cellfun,z)
  j = evaluate(gradient(phi),z)
  @notimplementedif D != Z
  # @fverdugo for the moment we assume positive jacobian
  # and that D == Z. To fixed soon.
  cellsum( f*(det(j)*w), dim=N )
end

function integrate(cellfun::CellFunction{Point{D},1},mesh::IntegrationMesh{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(mesh)
  integrate(cellfun,phi,quad)
end

function integrate(fun::Function,mesh::IntegrationMesh{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(mesh)
  cellfun = compose(fun,phi)
  integrate(cellfun,phi,quad)
end

end # module IntegrationMeshes
