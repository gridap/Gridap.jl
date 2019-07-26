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

tags = [7,6]
btrian = BoundaryTriangulation(model,tags)
quad = CellQuadrature(btrian,order=2)

bbh = restrict(bh,btrian)
cm = inner(bbh,bbh)

ca = integrate(cm,btrian,quad)
@test isa(ca,CellArray{Float64,2})
_ = collect(ca)

buh = restrict(uh,btrian)
cm = inner(bbh,buh)
ca = integrate(cm,btrian,quad)
@test isa(ca,CellArray{Float64,1})
_ = collect(ca)

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
order = 2
diritag = "boundary"
T = VectorValue{2,Float64}
fespace = ConformingFESpace(T,model,order,diritag)
ufun(x) = VectorValue(x[2],x[1])
uh = interpolate(fespace,ufun)
bh = FEBasis(fespace)
@test isa(ε(bh),FEBasis)
@test isa(inner(ε(bh),ε(bh)),CellMap{Point{2},1,Float64,3})
@test isa(inner(ε(bh),ε(uh)),CellMap{Point{2},1,Float64,2})

end # module
