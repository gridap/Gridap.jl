module CellBasesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.FESpaces

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

u = GenericCellBasis(Val{true}(),get_array(cell_basis),get_cell_map(cell_basis))

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

v = GenericCellBasis(Val{false}(),get_array(cell_basis),get_cell_map(cell_basis))

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

end # module
