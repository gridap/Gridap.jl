module FEFunctionsTests

using Gridap

using Test
using Gridap

using Gridap.FEFunctions: IndexCellFieldWithFEFunction

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

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

cf = IndexCellFieldWithFEFunction(uh.cellfield,uh)

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)
q = coordinates(quad)
v = collect(evaluate(cf,q))
g = collect(evaluate(gradient(cf),q))

@test isa(∇(cf),IndexCellFieldWithFEFunction)

test_index_cell_field(cf,q,v,g)

@test isa(∇(uh),IndexCellFieldWithFEFunction)
@test isa(Triangulation(∇(uh)),Triangulation)

u = CellField(trian,ufun)

e = u - uh
@test isa(e,IndexCellFieldWithFEFunction)

end # module
