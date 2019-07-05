module MultiFESpacesTests

using Test
using Gridap

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

ufun1(x) = x[1] + x[2]
U1 = TrialFESpace(fespace,ufun1)

ufun2(x) = x[1] + x[2]
U2 = TrialFESpace(fespace,ufun2)

U = MultiFESpace([U1,U2])

@test length(U) == 2
@test U[1] === U1
@test U[2] === U2

V1, state = iterate(U)
@test V1 === U1

V2, state = iterate(U,state)
@test V2 === U2

@test num_free_dofs(U) == num_free_dofs(U1) + num_free_dofs(U2)

a(v,u) = varinner(v,u)

bfun(x) = x[2] 
b(v) = varinner(v,CellField(trian,bfun))

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

v1 = CellBasis(V1)
v2 = CellBasis(V2)

u1 = CellBasis(U1)
u2 = CellBasis(U2)

mat11 = integrate(a(v1,u1),trian,quad)
mat12 = integrate(a(v1,u2),trian,quad)
mat22 = integrate(a(v2,u2),trian,quad)

mat = MultiCellMatrix([mat11,mat12,mat22],[(1,1),(1,2),(2,2)])

vec1 = integrate(b(v1),trian,quad)
vec2 = integrate(b(v2),trian,quad)

vec = MultiCellVector([vec1,vec2],[(1,),(2,)])

cellids = IdentityCellNumber(Int,length(vec1))
_vec = apply_constraints(U,vec,cellids)
@test _vec.blocks[1] === vec1
@test _vec.blocks[2] === vec2

_mat = apply_constraints_rows(U,mat,cellids)
@test _mat.blocks[1] === mat11
@test _mat.blocks[2] === mat12
@test _mat.blocks[3] === mat22

_mat = apply_constraints_cols(U,mat,cellids)
@test _mat.blocks[1] === mat11
@test _mat.blocks[2] === mat12
@test _mat.blocks[3] === mat22

dofs = celldofids(U)
@test dofs.blocks[1] === celldofids(U1)
@test dofs.blocks[2] === celldofids(U2)

n = num_free_dofs(U)
n1 = num_free_dofs(U[1])
n2 = num_free_dofs(U[2])
v = rand(n)

f1 = restrict_to_field(U,v,1)
@test f1 == v[1:n1]

f2 = restrict_to_field(U,v,2)
@test f2 == v[1+n1:n]

end # module MultiFESpacesTests
