module MultiFieldCellBasesTests

using Test
using Gridap.Geometry
using Gridap.Fields
using Gridap.Integration
using Gridap.FESpaces

using Gridap.MultiField
using Gridap.MultiField: CellBasisWithFieldID
using Gridap.MultiField: CellMatrixFieldWithFieldIds
using Gridap.MultiField: MultiFieldCellArray

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
order = 1

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)
q = get_coordinates(quad)

V = TestFESpace(model=model,reffe=:Lagrangian,order=order,conformity=:H1,valuetype=Float64)
Vp = TestFESpace(model=model,reffe=:Lagrangian,order=order,conformity=:H1,valuetype=Float64,dof_space=:physical)
U = TrialFESpace(V)

cell_basis_v = get_cell_basis(V)
field_id_v = 3
v = CellBasisWithFieldID(cell_basis_v,field_id_v)
vq = collect(evaluate(v,q))
test_cell_basis(v,q,vq)
@test is_in_ref_space(v)

field_id_a = 2
a = CellBasisWithFieldID(cell_basis_v,field_id_a)
aq = collect(evaluate(a,q))

cell_basis_u = get_cell_basis(U)
field_id_u = 4
u = CellBasisWithFieldID(cell_basis_u,field_id_u)
uq = collect(evaluate(u,q))

field_id_b = 5
b = CellBasisWithFieldID(cell_basis_u,field_id_b)
bq = collect(evaluate(b,q))

cell_basis_vp = get_cell_basis(Vp)
field_id_vp = 6
vp = CellBasisWithFieldID(cell_basis_vp,field_id_vp)
vpq = collect(evaluate(vp,q))
test_cell_basis(vp,q,vpq)
@test !is_in_ref_space(vp)

r = 2*v
test_cell_basis(r,q,2*vq)
@test isa(r,CellBasisWithFieldID)
@test r.field_id == field_id_v
@test is_test(r)

r = u*2
test_cell_basis(r,q,2*uq)
@test isa(r,CellBasisWithFieldID)
@test r.field_id == field_id_u
@test is_trial(r)

r = u + u
test_cell_basis(r,q,2*uq)
@test isa(r,CellBasisWithFieldID)
@test r.field_id == field_id_u
@test is_trial(r)

r = ∇(u)
@test isa(r,CellBasisWithFieldID)
@test r.field_id == field_id_u
@test is_trial(r)

r = v * u
test_cell_matrix_field(r,q,collect(evaluate(r,q)))
@test isa(r,CellMatrixFieldWithFieldIds)
@test r.field_id_rows == v.field_id
@test r.field_id_cols == u.field_id

r = u * v
test_cell_matrix_field(r,q,collect(evaluate(r,q)))
@test isa(r,CellMatrixFieldWithFieldIds)
@test r.field_id_rows == v.field_id
@test r.field_id_cols == u.field_id

s = a+v
@test s.blocks[1] === a
@test s.blocks[2] === v
@test s.block_ids[1] == (field_id_a,)
@test s.block_ids[2] == (field_id_v,)

s = a-v
@test s.blocks[1] === a
@test s.block_ids[1] == (field_id_a,)
@test s.block_ids[2] == (field_id_v,)

r = u*v
s = a*b

z = r + s
@test z.blocks[1] === r
@test z.blocks[2] === s
@test z.block_ids[1] == (field_id_v,field_id_u)
@test z.block_ids[2] == (field_id_a,field_id_b)

# Restrictions

btrian = BoundaryTriangulation(model)
bdegree = order
bquad = CellQuadrature(btrian,bdegree)

u_Γb = restrict(u,btrian)
@test isa(u_Γb,CellBasisWithFieldID)
@test is_trial(u_Γb)
@test u_Γb.field_id == field_id_u

strian = SkeletonTriangulation(model)
sdegree = order
squad = CellQuadrature(strian,sdegree)

u_Γs = restrict(u,strian)
@test isa(u_Γs.left,CellBasisWithFieldID)
@test isa(u_Γs.right,CellBasisWithFieldID)
@test u_Γs.left.field_id == field_id_u
@test u_Γs.right.field_id == field_id_u

jump_u_Γs = jump(u_Γs)
@test isa(jump_u_Γs.left,CellBasisWithFieldID)
@test isa(jump_u_Γs.right,CellBasisWithFieldID)
@test jump_u_Γs.left.field_id == field_id_u
@test jump_u_Γs.right.field_id == field_id_u

v_Γs = restrict(v,strian)

r = jump(u_Γs)*jump(v_Γs)
@test r.ll.field_id_rows == field_id_v
@test r.ll.field_id_cols == field_id_u

a_Γs = restrict(a,strian)
s = jump(u_Γs)*jump(a_Γs)
@test s.ll.field_id_rows == field_id_a
@test s.ll.field_id_cols == field_id_u

z = r + s
@test z.ll.block_ids == [(field_id_v,field_id_u), (field_id_a,field_id_u)]

# Integration

vec = integrate(v,trian,quad)
@test isa(vec,MultiFieldCellArray{Float64,1})

vec = integrate(v + a,trian,quad)
@test isa(vec,MultiFieldCellArray{Float64,1})

mat = integrate(v*u,trian,quad)
@test isa(mat,MultiFieldCellArray{Float64,2})

mat = integrate(v*u + a*b,trian,quad)
@test isa(mat,MultiFieldCellArray{Float64,2})

vec = integrate(jump(v_Γs),strian,squad)
@test isa(vec.left,MultiFieldCellArray{Float64,1})

vec = integrate(jump(v_Γs) + jump(a_Γs),strian,squad)
@test isa(vec.left,MultiFieldCellArray{Float64,1})
@test vec.left.block_ids == [(field_id_v,), (field_id_a,)]

mat = integrate(jump(v_Γs)*jump(u_Γs),strian,squad)
@test isa(mat.ll,MultiFieldCellArray{Float64,2})

mat = integrate(jump(v_Γs)*jump(u_Γs) + jump(a_Γs)*jump(u_Γs),strian,squad)
@test isa(mat.ll,MultiFieldCellArray{Float64,2})

end # module
