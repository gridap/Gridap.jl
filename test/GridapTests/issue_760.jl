module Issue760

using Gridap
using Test

cells = (2,2)
domain = (0,1,0,1)
model = CartesianDiscreteModel(domain,cells)

Ω = Interior(model)
Ω1 = Interior(model,[1,2])
Ω2 = Interior(model,[2,4])

f = CellField(1,Ω1)
x = get_cell_points(Ω2)
arr = f(x)
@test arr.maps.value.values_pos.value == [1,1,1,1]
@test arr.maps.value.values_neg.value == [0,0,0,0]

end # module
