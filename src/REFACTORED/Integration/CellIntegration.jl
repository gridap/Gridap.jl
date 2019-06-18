module CellIntegration

using Gridap

export integrate

function integrate(
  cellfun::CellMap{Point{D},1,T,N},
  phi::CellGeomap{D,Z},
  quad::CellQuadrature{D}) where {D,Z,T,N}
  z = coordinates(quad)
  w = weights(quad)
  f = evaluate(cellfun,z)
  j = evaluate(gradient(phi),z)
  # TODO this can be optimized wiht a Kernel
  cellsum( f*(meas(j)*w), dim=N )
end

function integrate(
  cellfun::CellMap{Point{D},1},trian::Triangulation{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = CellGeomap(trian)
  integrate(cellfun,phi,quad)
end

function integrate(
  fun::Function,trian::Triangulation{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = CellGeomap(trian)
  cellfun = compose(fun,phi)
  integrate(cellfun,phi,quad)
end

end # module CellIntegration
