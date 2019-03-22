
xe = Numa.DummyCellCoordinates2D(partition=(2,2))

basis = ShapeFunctionsScalarQua4()
cellbasis = Numa.CellBasisFromSingleInterpolation(basis)

quad = TensorProductQuadrature{2}(orders=[2,2])
cellquad = ConstantCellQuadrature(quad,length(xe))

q = coordinates(cellquad)

@test isa(cellbasis,Numa.EvaluableCellArray{2,Float64})

Ke = inner(cellbasis,cellbasis)

Keq = evaluate(Ke,q)

@test isa(Keq,CellMatrices{Float64})

gradcellbasis = gradient(cellbasis)

Ke = inner(gradcellbasis,gradcellbasis)

Keq = evaluate(Ke,q)

@test isa(Keq,CellMatrices{Float64})

