module CellIntegration

using Gridap
using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Polynomials
using Gridap.CellValues
using Gridap.CellMaps
using Gridap.CellQuadratures
using Gridap.Quadratures
using Gridap.Geometry

export integrate
import Gridap.CellMaps: CellField

function integrate(cellfun::CellMap{Point{D},1,T,N},phi::CellGeomap{D,Z},quad::CellQuadrature{D}) where {D,Z,T,N}
  z = coordinates(quad)
  w = weights(quad)
  f = evaluate(cellfun,z)
  j = evaluate(gradient(phi),z)
  cellsum( f*(meas(j)*w), dim=N )
end

function integrate(cellfun::CellMap{Point{D},1},trian::Triangulation{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(trian)
  integrate(cellfun,phi,quad)
end

function integrate(fun::Function,trian::Triangulation{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(trian)
  cellfun = compose(fun,phi)
  integrate(cellfun,phi,quad)
end

function CellField(trian::Triangulation,fun::Function)
  phi = geomap(trian)
  compose(fun,phi)
end

function CellField(trian::Triangulation{D,Z},fun::Function,u::CellField{Z}) where {D,Z}
  phi = geomap(trian)
  compose(fun,phi,u)
end

end # module CellIntegration
