module CellFieldsNewTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Integration
using Gridap.CellData
using Gridap.Fields
using Gridap.TensorValues
using LinearAlgebra
using FillArrays
using Test


node_to_x = 2*Point{2,Float64}[
  Point(0,0),
  Point(1,0),
  Point(2,0),
  Point(0,1),
  Point(1,1),
  Point(2,1),
  Point(0,2),
  Point(1,2),
  Point(2,2)]

cell_to_node =[
  [1,2,4,5],
  [2,3,5,6],
  [4,5,7,8],
  [5,6,8,9]]

cell_to_x = [ node_to_x[nodes] for nodes in cell_to_node]

l = length(cell_to_node)

reffe = QUAD4

degree = 2
refqua = Quadrature(get_polytope(reffe),degree)
q = Fill(get_coordinates(refqua),l)
w = Fill(get_weights(refqua),l)
npoin = num_points(refqua)

# cell map
ndofs = num_dofs(reffe)
s = GenericCellField( Fill(get_shapefuns(reffe),l), Fill((Base.OneTo(ndofs),),l), Val((:,)))
ϕ = lincomb(s,cell_to_x)
@test get_metasize(ϕ) == ()

@test ∇(ϕ) === ∇(ϕ)
@test get_metasize(∇(ϕ)) == ()
@test ∇(ϕ)(q) === ∇(ϕ)(q)
@test det(∇(ϕ)) === det(∇(ϕ))
@test det(∇(ϕ))(q) === det(∇(ϕ))(q)
r = apply(sum,det(∇(ϕ))(q))
test_array(r , Fill(4*npoin,l), ≈)

#In practice the inverse map of ϕ is never computed as long
# as you compose f∘ϕ before evaluating f
dv = s∘inverse_map(ϕ)
@test dv∘ϕ === s
@test (dv∘ϕ)(q) === (dv∘ϕ)(q)
@test ∇(dv) === ∇(dv)
@test ∇(dv)∘ϕ === ∇(dv)∘ϕ
@test (∇(dv)∘ϕ)(q) === (∇(dv)∘ϕ)(q)

# Gradient rules (test)
r1 = (∇(dv+dv)∘ϕ)(q)
r2 = ((∇(dv)+∇(dv))∘ϕ)(q)
test_array(r1,r2)

r1 = (∇(dv-dv)∘ϕ)(q)
r2 = ((∇(dv)-∇(dv))∘ϕ)(q)
test_array(r1,r2)

r1 = (∇(2*dv)∘ϕ)(q)
r2 = ((2*∇(dv))∘ϕ)(q)
test_array(r1,r2)

du = trialize_cell_basis(dv)

# Gradient rules (trial)
r1 = (∇(du+du)∘ϕ)(q)
r2 = ((∇(du)+∇(du))∘ϕ)(q)
test_array(r1,r2)

r1 = (∇(du-du)∘ϕ)(q)
r2 = ((∇(du)-∇(du))∘ϕ)(q)
test_array(r1,r2)

r1 = (∇(2*du)∘ϕ)(q)
r2 = ((2*∇(du))∘ϕ)(q)
test_array(r1,r2)

f(x) = x[2]
m = 3*∇(dv)⋅∇(du)*f + dv*du*2

dΩ_ref = CellQuadrature(q,w)
r = integrate(m,ϕ,dΩ_ref)
test_array(r,collect(r))

dΩ = ϕ(dΩ_ref)

r = ∫( 3*∇(dv)⋅∇(du)*f + dv*du*2 )*dΩ
test_array(r,collect(r))

#r = apply( Fields.IntKernel(), (m∘ϕ)(q), w, ∇(ϕ)(q))
#test_array(r,collect(r))

end # module
