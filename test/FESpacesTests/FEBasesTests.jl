module FEBasesTests

using Test
using Gridap

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

bh = FEBasis(fespace)

ufun(x) = x[1]

uh = interpolate(fespace,ufun)

@test isa(bh,FEBasis)
@test isa(+bh,FEBasis)
@test isa(-bh,FEBasis)
@test isa(∇(bh),FEBasis)
@test isa(bh+uh,FEBasis)
@test isa(uh-bh,FEBasis)
@test isa(inner(bh,bh),CellMap{Point{2},1,Float64,3})
@test isa(inner(bh,uh),CellMap{Point{2},1,Float64,2})

end # module
