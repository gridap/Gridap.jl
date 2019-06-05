module MultiFEBasesTests

using Test
using Gridap
using Gridap.FieldValues
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.CellMaps
using Gridap.Geometry.Cartesian
using Gridap.CellQuadratures
using Gridap.CellIntegration
using Gridap.MultiCellArrays
using Gridap.MultiFESpaces
using Gridap.MultiFEFunctions

import Gridap: gradient

u1fun(x) = x[1] + x[2]
u2fun(x) = x[1] - x[2]

u1fun_grad(x) = VectorValue(1.0,1.0)
u2fun_grad(x) = VectorValue(1.0,-1.0)

gradient(::typeof(u1fun)) = u1fun_grad
gradient(::typeof(u2fun)) = u2fun_grad

b1fun(x) = u2fun(x)
b2fun(x) = 0.0

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

b1field = CellField(trian,b1fun)
b2field = CellField(trian,b2fun)

order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

U1 = TrialFESpace(fespace,u1fun)

U2 = TrialFESpace(fespace,u2fun)

V1 = TestFESpace(fespace)
V2 = V1

U = MultiFESpace([U1,U2])
V = MultiFESpace([V1,V2])

u = FEBasis(U)
v = FEBasis(V)

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

end # module MultiFEBasesTests
