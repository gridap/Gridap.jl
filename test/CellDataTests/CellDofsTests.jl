module CellDofsTests

using Test
using FillArrays
using Gridap.Helpers
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Integration
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Geometry
using LinearAlgebra

domain = (0,1,0,1,0,1)
cells = (2,2,2)
model = simplexify(CartesianDiscreteModel(domain,cells))

trian = Triangulation(model)

v = GenericCellField(get_cell_shapefuns(trian),trian,ReferenceDomain())

ctype_to_reffe, cell_to_ctype = compress_cell_data(get_cell_reffe(trian))
ctype_to_dofbasis = map(get_dof_basis,ctype_to_reffe)
cell_dofbasis = expand_cell_data(ctype_to_dofbasis,cell_to_ctype)

s = CellDof(cell_dofbasis,trian,ReferenceDomain())

a = s(v)
r = fill(Matrix(I,4,4),num_cells(trian))
test_array(a,r)

end # module
