module MultiFieldCellFieldsTests

using Test
using Gridap.Geometry
using Gridap.CellData
using Gridap.MultiField

domain = (0,1,0,1)
cells = (2,2)
model = CartesianDiscreteModel(domain,cells)

trian = Triangulation(model)

u1(x) = sin(x[1])
cf1 = CellField(u1,trian)

u2(x) = cos(x[2])
cf2 = CellField(u2,trian)

cf = MultiFieldCellField([cf1,cf2])

@test cf1 === cf[1]
@test cf2 === cf[2]

_cf1, _cf2 = cf

@test cf1 === _cf1
@test cf2 === _cf2

end # module
