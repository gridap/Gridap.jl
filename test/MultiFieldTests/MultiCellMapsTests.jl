module MultiCellMapsTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.FieldValues
using Gridap.CellMaps
using Gridap.Geometry.Cartesian
using Gridap.FESpaces
using Gridap.CellIntegration
using Gridap.CellQuadratures
using Gridap.MultiCellArrays
using Gridap.MultiCellMaps
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

order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

V1 = TestFESpace(fespace)
V2 = V1

V = [V1,V2]
U = V

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

b1field = CellField(trian,b1fun)
b2field = CellField(trian,b2fun)

function a(v,u)
  a11 = varinner(∇(v[1]), ∇(u[1]))
  a12 = varinner(v[1],u[2])
  a22 = varinner(∇(v[2]), ∇(u[2]))
  return CellMap[a11,a12,a22], [(1,1),(1,2),(2,2)]
end

function b(v)
  b1 = varinner(v[1],b1field)
  b2 = varinner(v[2],b2field)
  return CellMap[b1,b2], [(1,),(2,)]
end

v = [ CellBasis(Vi) for Vi in V ]
u = [ CellBasis(Ui) for Ui in U ]

blocks, fieldids = a(v,u)

mcm = MultiCellMap(blocks,fieldids)

mca = integrate(mcm,trian,quad)

@test isa(mca,MultiCellArray{Float64,2})

u1 = CellField(trian,u1fun)
u2 = CellField(trian,u2fun)
u = [u1,u2]

blocks, fieldids = a(v,u)

fieldids = [ (fi[1],) for fi in fieldids ]

mcm = MultiCellMap(blocks,fieldids)
mcm1 = mcm

mca = integrate(mcm,trian,quad)

@test isa(mca,MultiCellArray{Float64,1})

blocks, fieldids = b(v)

mcm = MultiCellMap(blocks,fieldids)
mcm2 = mcm

mca = integrate(mcm,trian,quad)

@test isa(mca,MultiCellArray{Float64,1})

mcm = mcm1 + mcm2
mca = integrate(mcm,trian,quad)
@test isa(mca,MultiCellArray{Float64,1})

mcm = mcm1 - mcm2
mca = integrate(mcm,trian,quad)
@test isa(mca,MultiCellArray{Float64,1})

end
