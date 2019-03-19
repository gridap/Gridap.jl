
quad = TensorProductQuadrature{2}(orders=[2,2])
l = 10
cellquad = ConstantCellQuadrature(quad,l)

cellpoints = coordinates(cellquad)

basis = ShapeFunctionsScalarQua4()
cellbasis = Numa.CellBasisFromSingleInterpolation(basis)

cellbasisvalues = evaluate(cellbasis,cellpoints)

@test isa(cellbasis,CellBasis{Float64})

@test isa(cellbasisvalues,CellBasisValues{Float64})

@test length(cellbasisvalues) == l

cellbasisvalues2 = Numa.CellBasisValuesFromSingleInterpolation(basis,cellpoints)

@test isa(cellbasisvalues2,CellBasisValues{Float64})

for (values,values2) in zip(cellbasisvalues,cellbasisvalues2)
  @test values == values2
end

gradcellbasis = gradient(cellbasis)

@test isa(gradcellbasis,CellBasis{VectorValue{2}})

gradcellbasisvalues = evaluate(gradcellbasis,cellpoints)

@test isa(gradcellbasisvalues,CellBasisValues{VectorValue{2}})

