
imesh = Numa.DummyIntegrationMesh2D(partition=(3,3))

@test isa(imesh,Numa.IntegrationMesh)

coords = cellcoordinates(imesh)

basis = Numa.cellbasis(imesh)

phi = geomap(imesh)

@test isa(coords,CellPoints{2})

@test isa(basis,CellBasis{2,Float64})

@test isa(phi,CellField{2,Point{2}})

