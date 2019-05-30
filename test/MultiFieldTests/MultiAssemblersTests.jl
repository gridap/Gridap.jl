include("../../src/MultiField/MultiAssemblers.jl")

module MultiAssemblersTests

using Test

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Assemblers
using Gridap.Geometry.Cartesian
using Gridap.CellMaps
using Gridap.MultiCellArrays
using Gridap.CellQuadratures
using Gridap.CellIntegration

using ..MultiAssemblers

order = 1
diritag = "boundary"
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(2,3))
fespace = ConformingFESpace(Float64,model,order,diritag)

ufun1(x) = x[1] + x[2]
U1 = TrialFESpace(fespace,ufun1)

ufun2(x) = x[1] + x[2]
U2 = TrialFESpace(fespace,ufun2)

V1 = TestFESpace(fespace)
V2 = V1

V = [V1,V2]
U = [U1,U2]

assem = SparseMatrixAssembler(V,U)

@test isa(assem,MultiSparseMatrixAssembler)

a(v,u) = varinner(∇(v),∇(u))

m(v,u) = varinner(v,u)

bfun1(x) = 1.0
b1(v) = varinner(v,CellField(trian,bfun1))

bfun2(x) = 2.0
b2(v) = varinner(v,CellField(trian,bfun2))

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

v1 = CellBasis(V1)
v2 = CellBasis(V2)

u1 = CellBasis(U1)
u2 = CellBasis(U2)

mat11 = integrate(a(v1,u1),trian,quad)
mat12 = integrate(m(v1,u2),trian,quad)
mat22 = integrate(m(v2,u2),trian,quad)

mat = MultiCellMatrix([mat11,mat12,mat22],[(1,1),(1,2),(2,2)])

vec1 = integrate(b1(v1),trian,quad)
vec2 = integrate(b2(v2),trian,quad)

vec = MultiCellVector([vec1,vec2],[(1,),(2,)])

v = assemble(assem,vec)
@test v ≈ [0.166666666666, 0.166666666666, 0.3333333333333, 0.333333333333]

mm = assemble(assem,mat)

assemble!(mm,assem,mat)

end # module MultiAssemblersTests
