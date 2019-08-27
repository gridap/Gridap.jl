module FEFunctionsTests

using Gridap

using Test
using Gridap

using Gridap.Triangulations: IndexCellFieldWithTriangulation

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

trian = Triangulation(model)

ufun(x) = x[1]

uh = interpolate(fespace,ufun)

@test isa(Triangulation(uh),Triangulation)

tags = [7,6]
btrian = BoundaryTriangulation(model,tags)
bquad = CellQuadrature(btrian,order=0)

buh = restrict(uh,btrian)

cn = integrate(buh,btrian,bquad)
@test isa(cn,CellNumber)
_ = collect(cn)

cf = IndexCellFieldWithTriangulation(uh.cellfield,trian)

quad = CellQuadrature(trian,order=2)
q = coordinates(quad)
v = collect(evaluate(cf,q))
g = collect(evaluate(gradient(cf),q))

@test isa(∇(cf),IndexCellFieldWithTriangulation)

test_index_cell_field(cf,q,v,g)

@test isa(∇(uh),IndexCellFieldWithTriangulation)
@test isa(Triangulation(∇(uh)),Triangulation)

u = CellField(trian,ufun)

@test isa(u,IndexCellFieldWithTriangulation)

e = u - uh
@test isa(e,IndexCellFieldWithTriangulation)

end # module
