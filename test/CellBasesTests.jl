
imesh = Numa.DummyIntegrationMesh2D(partition=(2,2))

quad = TensorProductQuadrature{2}(orders=[2,2])
l = ncells(imesh)
cellquad = ConstantCellQuadrature(quad,l)

cellpoints = coordinates(cellquad)

basis = ShapeFunctionsScalarQua4()
cellbasis = Numa.CellBasisFromSingleInterpolation(basis)

cellbasisvalues = evaluate(cellbasis,cellpoints)

@test isa(cellbasis,CellBasis{2,Float64})

@test isa(cellbasisvalues,CellBasisValues{Float64})

@test length(cellbasisvalues) == l

cellbasisvalues2 = Numa.CellBasisValuesFromSingleInterpolation(basis,cellpoints)

@test isa(cellbasisvalues2,CellBasisValues{Float64})

for (values,values2) in zip(cellbasisvalues,cellbasisvalues2)
  @test values == values2
end

gradcellbasis = gradient(cellbasis)

@test isa(gradcellbasis,CellBasis{2,VectorValue{2}})

gradcellbasisvalues = evaluate(gradcellbasis,cellpoints)

@test isa(gradcellbasisvalues,CellBasisValues{VectorValue{2}})

phi = geomap(imesh)

physcellbasis = mapderivatives(cellbasis,phi)

physgradcellbasis = gradient(physcellbasis)

physgradcellbasis_q = evaluate(physgradcellbasis,cellpoints)

for dNcg in physgradcellbasis_q
  for dNg in dNcg
    @test isa(dNg,VectorValue{2})
  end
end

