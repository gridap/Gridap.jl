module FEBasesTests

using Test
using Gridap

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

bh = FEBasis(fespace)

@test isa(Triangulation(bh),Triangulation)

ufun(x) = x[1]

uh = interpolate(fespace,ufun)

@test isa(bh,FEBasis)
@test isa(+bh,FEBasis)
@test isa(-bh,FEBasis)
@test isa(∇(bh),FEBasis)
@test isa(bh+uh,FEBasis)
@test isa(uh-bh,FEBasis)
@test isa(bh-1.0,FEBasis)
@test isa(1.0+bh,FEBasis)
@test isa(inner(bh,bh),CellMap{Point{2},1,Float64,3})
@test isa(inner(bh,uh),CellMap{Point{2},1,Float64,2})

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

σfun(x,u,i) = x[1] + u
σ(u,i) = CellBasis(trian,σfun,u,i)
ids = ones(Int,ncells(trian))
cm = inner(bh,σ(bh,ids))
m = integrate(cm,trian,quad)
_ = collect(m)

tags = [7,6]
btrian = BoundaryTriangulation(model,tags)
bquad = CellQuadrature(btrian,order=2)

bbh = restrict(bh,btrian)
cm = inner(bbh,bbh)

ca = integrate(cm,btrian,bquad)
@test isa(ca,CellArray{Float64,2})
_ = collect(ca)

buh = restrict(uh,btrian)
cm = inner(bbh,buh)
ca = integrate(cm,btrian,bquad)
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
@test isa(div(bh),FEBasis)
@test isa(curl(bh),FEBasis)
@test isa(inner(ε(bh),ε(bh)),CellMap{Point{2},1,Float64,3})
@test isa(inner(ε(bh),ε(uh)),CellMap{Point{2},1,Float64,2})

σfun(x,ε) = 4*tr(ε)*one(ε) + 2*ε
σ(u) = CellBasis(trian,σfun,u)
cm = inner(ε(bh),σ(ε(bh)))
m = integrate(cm,trian,quad)
_ = collect(m)

end # module
