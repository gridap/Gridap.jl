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
ϕ = GenericCellMap(lincomb(fl,cl))

degree = 3
quad = CellQuadrature(degree,[QUAD,],Fill(1,ncells))
q = get_coordinates(quad)
dv = GenericCellField(Fill(b,ncells),Fill((Base.OneTo(ndofs),),ncells),Val((:,)))∘inverse_map(ϕ)
u = convert_to_cell_field(x->x[2]^2,length(ϕ))

@law g(u) = u^2
@law h(u,dv) = (1+u)*dv

cf = g(u)
@test isa(cf,CellField)
@test get_metasize(cf) == ()
test_array(evaluate(cf∘ϕ,q),evaluate(operate(g,u)∘ϕ,q),≈)

cf = h(u,dv)
@test isa(cf,CellField)
@test get_metasize(cf) == (:,)
test_array(evaluate(cf∘ϕ,q),evaluate(operate(h,u,dv)∘ϕ,q),≈)


end # module
