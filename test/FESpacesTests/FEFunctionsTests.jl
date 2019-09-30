module FEFunctionsTests

using Gridap

using Test
using Gridap

using Gridap.Triangulations: IndexCellFieldWithTriangulation

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
order = 1
diritag = "boundary"
fespace = H1ConformingFESpace(Float64,model,order,diritag)

trian = Triangulation(model)

ufun(x) = x[1]

uh = interpolate(fespace,ufun)

@test string(uh) == "FEFunction object"
sr = "FEFunction object:\n physdim: 2\n refdim: 2\n valuetype: Float64\n doftype: Float64\n nfree: 9\n ndiri: 16\n ncells: 16"
@test sprint(show,"text/plain",uh) == sr

@test isa(Triangulation(uh),Triangulation)

tags = [7,6]
btrian = BoundaryTriangulation(model,tags)
bquad = CellQuadrature(btrian,degree=0)

buh = restrict(uh,btrian)

cn = integrate(buh,btrian,bquad)
@test isa(cn,CellNumber)
_ = collect(cn)

cf = IndexCellFieldWithTriangulation(uh.cellfield,trian)

quad = CellQuadrature(trian,degree=2)
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

e = ufun - uh
@test isa(e,IndexCellFieldWithTriangulation)

e = uh - ufun
@test isa(e,IndexCellFieldWithTriangulation)

T = VectorValue{2,Float64}
fespace = H1ConformingFESpace(T,model,order,diritag)
ufun(x) = VectorValue(x[2],x[1])
uh = interpolate(fespace,ufun)

@test isa(div(uh),IndexCellFieldWithTriangulation)
@test isa(curl(uh),IndexCellFieldWithTriangulation)



end # module
