module MultiFEBasesTests

using Test
using Gridap
using Gridap.MultiFEBases: FEBasisWithFieldId

import Gridap: ∇

u1fun(x) = x[1] + x[2]
u2fun(x) = x[1] - x[2]

u1fun_grad(x) = VectorValue(1.0,1.0)
u2fun_grad(x) = VectorValue(1.0,-1.0)

∇(::typeof(u1fun)) = u1fun_grad
∇(::typeof(u2fun)) = u2fun_grad

b1fun(x) = u2fun(x)
b2fun(x) = 0.0

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

trian = Triangulation(model)
quad = CellQuadrature(trian,degree=2)

b1field = CellField(trian,b1fun)
b2field = CellField(trian,b2fun)

order = 1
diritag = "boundary"
fespace = H1ConformingFESpace(Float64,model,order,diritag)

U1 = TrialFESpace(fespace,u1fun)

U2 = TrialFESpace(fespace,u2fun)

V1 = TestFESpace(fespace)
V2 = V1

U = MultiFESpace([U1,U2])
V = MultiFESpace([V1,V2])

u = FEBasis(U)
v = FEBasis(V)

@test isa(u[1]+1,FEBasisWithFieldId)
@test isa(1-u[1],FEBasisWithFieldId)

@test isa(Triangulation(u[1]),Triangulation)

a(v,u) = inner(∇(v[1]),∇(u[1])) + inner(v[1],u[2]) + inner(∇(v[2]),∇(u[2]))

b(v) = inner(v[1],b1field) + inner(v[2],b2field)

mcm = a(v,u)

mca = integrate(mcm,trian,quad)

@test isa(mca,MultiCellArray{Float64,2})

mcm = b(v)

mca = integrate(mcm,trian,quad)

@test isa(mca,MultiCellArray{Float64,1})

u1 = CellField(trian,u1fun)
u2 = CellField(trian,u2fun)
u = [u1,u2]

mcm = a(v,u)

mca = integrate(mcm,trian,quad)

@test isa(mca,MultiCellArray{Float64,1})

zh = zero(U)

mcm = a(v,zh)

mca = integrate(mcm,trian,quad)

@test isa(mca,MultiCellArray{Float64,1})

T = VectorValue{2,Float64}
order = 1
diritag = "boundary"
fespace = H1ConformingFESpace(T,model,order,diritag)

V1 = TestFESpace(fespace)
V2 = V1

V = MultiFESpace([V1,V2])

bh = FEBasis(V)

@test isa(∇(bh[1]),FEBasisWithFieldId)
@test isa(ε(bh[1]),FEBasisWithFieldId)
@test isa(div(bh[1]),FEBasisWithFieldId)
@test isa(curl(bh[1]),FEBasisWithFieldId)

σfun(x,u,i) = i*u
ids = ones(Int,ncells(trian))
cb = CellBasis(trian,σfun,bh[1],ids)
@test isa(cb,FEBasisWithFieldId)

strian = SkeletonTriangulation(model)
squad = CellQuadrature(strian,degree=2)

sbh = restrict(bh,strian)
@test isa(sbh[1].cellfield1,FEBasisWithFieldId)
@test isa(sbh[1].cellfield2,FEBasisWithFieldId)

cm = inner(jump(sbh[1]),jump(sbh[1]))
ca = integrate(cm,strian,squad)

@test isa(ca.cellmatrix11,MultiCellArray)

uh = zero(V)
suh = restrict(uh,strian)
@test isa(suh[1].cellfield1,CellField)
@test isa(suh[1].cellfield2,CellField)

cm = inner(jump(sbh[1]),jump(suh[1]))
ca = integrate(cm,strian,squad)

end # module MultiFEBasesTests
