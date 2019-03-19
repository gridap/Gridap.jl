
abstract type IntegrationMesh{D} end

cellcoordinates(::IntegrationMesh{D} where D)::CellPoints{D} = @abstractmethod

cellbasis(::IntegrationMesh)::CellBasis{Float64} = @abstractmethod

# TODO
#function geomap(self::IntegrationMesh{D}) where D
#  coords = cellcoordinates(self)
#  basis = cellbasis(self)
#  interpolate(coords,basis)
#end
