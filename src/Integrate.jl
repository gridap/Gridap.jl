
function integrate(
  f::EvaluableCellArray{Z,T,N}, dom::IntegrationDomain{Z,D},
  quad::CellQuadrature{Z}) where {D,T,N,Z}

  phi = geomap(dom)
  jac = gradient(phi)
  pointsZ = coordinates(quad)
  weightsZ = weights(quad)
  fZ = evaluate(f,pointsZ)
  cellsum( fZ * weightsZ * det(jac), dim=N )

end

