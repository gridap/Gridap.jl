module CellFieldsTests

#This is a first steps to decouple CellFields from the Geometry

using Test
using BlockArrays
using FillArrays
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: IntMap
using Gridap.Fields: FieldOpMap
using Gridap.Fields: trialize_basis_value
using Gridap.TensorValues
using Gridap.Fields: MockField, MockBasis, OtherMockBasis

# In this furure all this shoul be exported from Fields
using Gridap.Geometry: GenericCellField, get_cell_axes
using Gridap.Geometry: convert_to_cell_field, get_metasize
using Gridap.Geometry: CellField, similar_object, trialize_cell_basis
using Gridap.Geometry: merge_skeleton_dof_ids, merge_skeleton_cell_fields
using Gridap.Geometry: is_basis, is_test, is_trial

np = 3
ndofs = 4

p = Point(1,2)
x = fill(p,np)
z = 2.0

v = VectorValue(3.0,1.5)
w = VectorValue(3.4,3.5)
a = MockBasis{2}(v,ndofs)
b = MockBasis{2}(w,ndofs)
c = fill(1.0,ndofs)
f = OtherMockBasis{2}(ndofs)

g = MockField{2}(v)

l = 10
xl = Fill(x,l)
zl = [ z for  i in 1:l]
cl = fill(c,l)
fl = Fill(f,l)
ϕl = lincomb(fl,cl)
gl = fill(g,l)
al = Fill(a,l)
bl = fill(b,l)

gf = GenericCellField(gl,ϕl,Val(true))
af = GenericCellField(al,ϕl,Val(true),Fill((Base.OneTo(ndofs),),l),Val((:,)))
bf = GenericCellField(bl,ϕl,Val(true),Fill((Base.OneTo(ndofs),),l),Val((:,)))
zf = convert_to_cell_field(zl,ϕl)
df = af*zf
dft = trialize_cell_basis(df)

# Check memoization
df_x1 = evaluate(df,xl)
df_x2 = evaluate(df,xl)
@test df_x1 === df_x2
∇gf1 = ∇(gf)
∇gf2 = ∇(gf)
@test ∇gf1 === ∇gf2
@test evaluate(∇gf1,xl) === evaluate(∇gf2,xl)
εgf1 = ε(gf)
εgf2 = ε(gf)
@test εgf1 === εgf2
@test ∇×gf === ∇×gf
@test evaluate(∇×gf,xl) === evaluate(∇×gf,xl)

@test is_test(af)
@test is_trial(dft)
@test is_basis(af)
@test is_basis(dft)
mf = af⋅dft
@test get_metasize(mf) == (:,:)
mf_x = evaluate(mf,xl)
@test size(mf_x[1]) == (np,ndofs,ndofs)

@test get_cell_axes(mf) == Fill((Base.OneTo(ndofs),Base.OneTo(ndofs)),l)

idsL = [ i*collect(1:ndofs)  for i in 1:l]
idsR = [ 2*i*collect(1:ndofs)  for i in 1:l]
axesL = Fill((Base.OneTo(ndofs),),l)
axesR = Fill((Base.OneTo(ndofs),),l)

idsS = merge_skeleton_dof_ids(idsL,idsR,axesL,axesR)
@test isa(idsS,VectorOfBlockArrayCoo)

afL, afR = merge_skeleton_cell_fields(af,2*af)

afL_x = evaluate(afL,xl)
afR_x = evaluate(afR,xl)
@test isa(afL_x,VectorOfBlockArrayCoo)
@test isa(afR_x,VectorOfBlockArrayCoo)

end # module
