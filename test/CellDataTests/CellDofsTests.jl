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

x = get_cell_points(s)
vx = v(x)
test_array(vx,r)

acell_to_cell = [1,3,5,3]
atrian = Triangulation(trian,acell_to_cell)
s_a = change_domain(s,atrian,DomainStyle(s))
v_a = change_domain(v,atrian,DomainStyle(v))

a = s_a(v_a)
r = fill(Matrix(I,4,4),num_cells(atrian))
test_array(a,r)

x_a = get_cell_points(s_a)
vx_a = v(x_a)
test_array(vx_a,r)

end # module
