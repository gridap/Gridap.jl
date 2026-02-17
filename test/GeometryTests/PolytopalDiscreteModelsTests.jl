module PolytopalDiscreteModelsTests
using Gridap
using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Fields, Gridap.FESpaces, Gridap.CellData

using FillArrays

function test_model(model)
  D = num_cell_dims(model)

  # Model API
  test_discrete_model(model)
  Geometry.restrict(model,[1,2,3,4])

  # Triangulations

  for d in 0:D
    trian = Triangulation(ReferenceFE{d},model)
    test_triangulation(trian)
  end

  Ω = Triangulation(ReferenceFE{D},model,[1,2,3,4])
  Γ = Boundary(model)
  Λ = Skeleton(model)

  d = tempdir()
  writevtk(model,d*"/polytopal_model")
  writevtk(Ω,d*"/polytopal_trian")
  writevtk(Γ,d*"/polytopal_boundary")
  writevtk(Λ,d*"/polytopal_skeleton")
end

model = CartesianDiscreteModel((0,1,0,1),(2,2))

pmodel = Geometry.PolytopalDiscreteModel(model)
pgrid = Geometry.PolytopalGrid(get_polytopes(pmodel))
test_model(pmodel)

vmodel = Geometry.voronoi(Geometry.simplexify(model))
test_model(vmodel)

model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
pmodel = Geometry.PolytopalDiscreteModel(model)
test_model(pmodel)


end