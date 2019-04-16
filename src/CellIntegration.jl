module CellIntegration

using Numa
using Numa.Helpers
using Numa.FieldValues
using Numa.Polynomials
using Numa.CellValues
using Numa.CellFunctions
using Numa.CellQuadratures
using Numa.Quadratures
using Numa.Geometry

export integrate
export cellfield

function integrate(cellfun::CellFunction{Point{D},1,T,N},phi::CellGeomap{D,Z},quad::CellQuadrature{D}) where {D,Z,T,N}
  z = coordinates(quad)
  w = weights(quad)
  f = evaluate(cellfun,z)
  j = evaluate(gradient(phi),z)
  cellsum( f*(meas(j)*w), dim=N )
end

function integrate(cellfun::CellFunction{Point{D},1},trian::Triangulation{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(trian)
  integrate(cellfun,phi,quad)
end

function integrate(fun::Function,trian::Triangulation{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(trian)
  cellfun = compose(fun,phi)
  integrate(cellfun,phi,quad)
end

function cellfield(trian::Triangulation,fun::Function)
  phi = geomap(trian)
  compose(fun,phi)
end

function cellfield(trian::Triangulation{D,Z},fun::Function,u::CellField{Z}) where {D,Z}
  phi = geomap(trian)
  compose(fun,phi,u)
end

end # module CellIntegration
