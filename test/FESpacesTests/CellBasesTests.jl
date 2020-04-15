module CellBasesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.FESpaces
using Gridap.FESpaces: SkeletonCellBasis
using Gridap.FESpaces: ReducedSkeletonCellBasis
using Gridap.FESpaces: SkeletonCellVector
using Gridap.FESpaces: SkeletonCellMatrixField
using Gridap.FESpaces: SkeletonCellMatrix

using Gridap.Fields: MockField, MockBasis

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

dirichlet_tags = ["tag_1","tag_6"]
V0 = GradConformingFESpace(reffes,model,dirichlet_tags)

trian = get_triangulation(model)

degree = 2
quad = CellQuadrature(trian,degree)
q = get_coordinates(quad)

cell_basis = get_cell_basis(V0)

cbx = collect(evaluate(cell_basis,q))

f(x) = sin(4*pi*(x[1]-x[2]^2))+1

uh = interpolate_everywhere(V0,f)
uhx = collect(evaluate(uh,q))

v = GenericCellBasis(Val{false}(),get_array(cell_basis),get_cell_map(cell_basis),Val{true}())
u = GenericCellBasis(Val{true}(),get_array(cell_basis),get_cell_map(cell_basis),Val{true}())

dv = v
du = u

a = rand(length(v))

r = operate(+,uh,2*uh,a)
@test isa(r,CellField)
rq = collect(evaluate(r,q))

r = operate(+,du,uh,a)
@test isa(r,CellBasis)
@test is_trial(r)
rq = collect(evaluate(r,q))

r = operate(+,uh,du,a)
@test isa(r,CellBasis)
@test is_trial(r)
rq = collect(evaluate(r,q))

r = operate(+,du,du,uh,a)
@test isa(r,CellBasis)
@test is_trial(r)
rq = collect(evaluate(r,q))

r = u*2
test_cell_basis(r,q,2*cbx)
@test is_trial(r)
@test ! is_test(r)

r = 2*u
test_cell_basis(r,q,2*cbx)
@test is_trial(r)

r = u + uh
rr = [ ai .+ bi for (ai,bi) in  zip(cbx,uhx)]
test_cell_basis(r,q,rr)
@test is_trial(r)

r = v*2
test_cell_basis(r,q,2*cbx)
@test is_test(r)
@test ! is_trial(r)

r = 2*v
test_cell_basis(r,q,2*cbx)
@test is_test(r)

r = v + uh
rr = [ ai .+ bi for (ai,bi) in  zip(cbx,uhx)]
test_cell_basis(r,q,rr)
@test is_test(r)

w = inner(u,v)
@test isa(w,CellMatrixField)

w = u * v
@test isa(w,CellMatrixField)
wx = collect(evaluate(w,q))

s = 2*w
test_cell_matrix_field(s,q,2*wx,≈)

s = w*2
test_cell_matrix_field(s,q,2*wx,≈)

s = w+w
test_cell_matrix_field(s,q,2*wx,≈)

s = (3*w)-w
test_cell_matrix_field(s,q,2*wx,≈)

mat = integrate(w,trian,quad)
mat = collect(mat)

vec = integrate(v,trian,quad)
vec = collect(vec)

mat2 = integrate(2*w,trian,quad)
test_array(mat2,2*mat)

vec2 = integrate(2*v,trian,quad)
test_array(vec2,2*vec)

# Testing boundary

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,degree)
bq = get_coordinates(bquad)

bv = restrict(v,btrian)
@test isa(bv,CellBasis)
@test is_test(bv)
bvec = integrate(bv,btrian,bquad)
bvec = collect(bvec)

bvec2 = integrate(2*bv,btrian,bquad)
test_array(bvec2,2*bvec)

bu = restrict(u,btrian)
@test isa(bv,CellBasis)
@test is_trial(bu)

bw = bv * bu
@test isa(bw,CellMatrixField)
bmat = integrate(bw,btrian,bquad)
bmat = collect(bmat)

bmat2 = integrate(bw*2,btrian,bquad)
test_array(bmat2,2*bmat)

# Testing Skeleton

strian = SkeletonTriangulation(model)
squad = CellQuadrature(strian,degree)
sq = get_coordinates(squad)

sv = restrict(v,strian)
@test isa(sv,SkeletonCellBasis)
@test is_test(sv)
@test sv.left === sv.inward
@test sv.right === sv.outward

_sv = jump(2*∇(sv))
@test isa(_sv,ReducedSkeletonCellBasis)
@test is_test(_sv)

_sv = jump(∇(sv)*f)
@test isa(_sv,ReducedSkeletonCellBasis)
@test is_test(_sv)

sv = jump(sv)
@test isa(sv,ReducedSkeletonCellBasis)
@test is_test(sv)

su = restrict(u,strian)
@test isa(su,SkeletonCellBasis)
@test is_trial(su)

su = jump(su)
@test isa(su,ReducedSkeletonCellBasis)
@test is_trial(su)

svec = integrate(sv,strian,squad)
@test isa(svec,SkeletonCellVector)
svec_left = collect(svec.left)
svec_right = collect(svec.right)

svec2 = integrate(2*sv,strian,squad)
test_array(svec2.left,2*svec_left)
test_array(svec2.right,2*svec_right)

sw = sv * su
@test isa(sw,SkeletonCellMatrixField)

smat = integrate(sw,strian,squad)
@test isa(smat,SkeletonCellMatrix)

smat_ll = collect(smat.ll)
smat_lr = collect(smat.lr)
smat_rl = collect(smat.rl)
smat_rr = collect(smat.rr)

smat2 = integrate(2*sw,strian,squad)
test_array(smat2.ll,2*smat_ll)
test_array(smat2.lr,2*smat_lr)
test_array(smat2.rl,2*smat_rl)
test_array(smat2.rr,2*smat_rr)

end # module
