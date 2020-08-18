module LawTests

using Test
using FillArrays
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.CellData
using Gridap.Fields: OtherMockBasis

ndofs = num_dofs(QUAD4)
b = get_shapefuns(QUAD4)
c = fill(1,ndofs)
f = OtherMockBasis{2}(ndofs)

ncells = 10
cl = fill(c,ncells)
fl = Fill(f,ncells)
ϕl = lincomb(fl,cl)

degree = 3
quad = CellQuadrature(degree,[QUAD,],Fill(1,ncells))
q = get_coordinates(quad)
dv = GenericCellField(Fill(b,ncells),ϕl,Val(true),Fill((Base.OneTo(ndofs),),ncells),Val((:,)))
u = convert_to_cell_field(x->x[2]^2,ϕl)

@law g(u) = u^2
@law h(u,dv) = (1+u)*dv

cf = g(u)
@test isa(cf,CellField)
@test get_metasize(cf) == ()
test_array(evaluate(cf,q),evaluate(operate(g,u),q),≈)

cf = h(u,dv)
@test isa(cf,CellField)
@test get_metasize(cf) == (:,)
test_array(evaluate(cf,q),evaluate(operate(h,u,dv),q),≈)


end # module
