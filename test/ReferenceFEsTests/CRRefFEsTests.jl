using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Helpers

using FillArrays


p = TRI
D = num_dims(p)

# T = Float64
T = VectorValue{D,Float64}


cr_reffe = CRRefFE(T,p,1)

cr_dofs = get_dof_basis(cr_reffe)
cr_prebasis  = get_prebasis(cr_reffe)
cr_shapefuns = get_shapefuns(cr_reffe) 

M = evaluate(cr_dofs, cr_shapefuns)

partition = (0,1,0,1)
cells = (2,2)
model = simplexify(CartesianDiscreteModel(partition, cells))

V = FESpace(model,cr_reffe)
get_cell_dof_ids(V)