module CellDofBasesTests

using Test
using FillArrays
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Fields: MockField, OtherMockBasis
using LinearAlgebra

reffe = QUAD4
ndofs = 3
b = get_shapefuns(reffe)
d = get_dof_basis(reffe)
c = fill(1.0,ndofs)
f = OtherMockBasis{2}(ndofs)

l = 10
bl = fill(b,l)
dl = Fill(d,l)
cl = fill(c,l)
fl = Fill(f,l)
ϕ = GenericCellMap(lincomb(fl,cl))

ref_style = Val(true)
axs = (Base.OneTo(num_dofs(reffe)),)
cell_axes = Fill(axs,l)
cb = GenericCellField(bl,cell_axes,Val((:,)))
cd = GenericCellDof(dl)
test_cell_dof_basis(cd,cb)

@test evaluate(cd,cb) ≈ Fill(Matrix(I,num_dofs(reffe),num_dofs(reffe)),l)

cb_phys = cb∘inverse_map(ϕ)
cd_phys = ϕ(cd)
@test cd_phys(cb_phys) ≈ Fill(Matrix(I,num_dofs(reffe),num_dofs(reffe)),l)

end # module
