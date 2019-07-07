module FEFunctionsTests

using Gridap

using Test
using Gridap

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

ufun(x) = x[1]

uh = interpolate(fespace,ufun)

tags = [7,6]
btrian = BoundaryTriangulation(model,tags)
quad = CellQuadrature(btrian,order=0)

buh = restrict(uh,btrian)

cn = integrate(buh,btrian,quad)
@test isa(cn,CellNumber)
_ = collect(cn)

end # module
