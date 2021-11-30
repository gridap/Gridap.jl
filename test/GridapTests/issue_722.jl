module Issue722

using Test
using Gridap
using Gridap.Fields
using Gridap.CellData

model = simplexify(CartesianDiscreteModel((0,1,0,1),(2,2)))
Ω=Triangulation(model)
Ωo=Gridap.Geometry.TriangulationView(Ω,[1,2,3])
glue = get_glue(Ωo,Val(2))
@test glue.tface_to_mface == [1, 2, 3]
@test glue.mface_to_tface == [1, 2, 3, -1, -2, -3, -4, -5]
no=GenericCellField(fill(GenericField(x->x[2]),num_cells(Ωo)),Ωo,ReferenceDomain())
x_Ω = get_cell_points(Ω)
nox = no(x_Ω)
collect(nox)

end # module
