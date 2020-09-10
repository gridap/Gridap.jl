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
q = GenericCellPoint(Fill(get_coordinates(refqua),l))
w = Fill(get_weights(refqua),l)
npoin = num_points(refqua)

# cell map
ndofs = num_dofs(reffe)
s = GenericCellField( Fill(get_shapefuns(reffe),l), Fill((Base.OneTo(ndofs),),l), Val((:,)))
ϕ = GenericCellMap(get_array(lincomb(s,cell_to_x)))
x = ϕ(q)

test_cell_point(q)
test_cell_map(ϕ)
test_cell_point(x)

@test isa(∇(ϕ),CellField)
@test ∇(ϕ) === ∇(ϕ)
@test get_metasize(∇(ϕ)) == ()
@test ϕ(q) === ϕ(q)
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

cell_to_u = [ rand(length(nodes)) for nodes in cell_to_node]
uh = lincomb(du,cell_to_u)
r = (uh∘ϕ)(q)
test_array(r,collect(r))

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

# Gradient rules (FEFunction)
r1 = (∇(uh)∘ϕ)(q)
r2 = (lincomb(∇(du),cell_to_u)∘ϕ)(q)
test_array(r1,r2)

r2 = (lincomb(∇(dv)∘ϕ,cell_to_u))(q)
test_array(r1,r2)

r1 = (∇(uh+uh)∘ϕ)(q)
r2 = ((∇(uh)+∇(uh))∘ϕ)(q)
test_array(r1,r2)

r1 = (∇(uh-uh)∘ϕ)(q)
r2 = ((∇(uh)-∇(uh))∘ϕ)(q)
test_array(r1,r2)

r1 = (∇(2*uh)∘ϕ)(q)
r2 = ((2*∇(uh))∘ϕ)(q)
test_array(r1,r2)

# integration

f(x) = x[2]
m = 3*∇(dv)⋅∇(du)*f + dv*du*2

dΩ_ref = CellQuadrature(q,w)
r = integrate(m,ϕ,dΩ_ref)
test_array(r,collect(r))

dΩ = ϕ(dΩ_ref)

r = ∫( 3*∇(dv)⋅∇(du)*f + dv*du*2 )*dΩ
test_array(r,collect(r))

# FaceMap

cell_map = ϕ
face_to_cell = [2,3,1]
cf = convert_to_cell_field(identity,length(face_to_cell))
refface_to_refcell_map = GenericCellMap(get_array(cf),face_to_cell,l)

ϕ_Γ_1 = GenericCellMap(reindex(get_array(ϕ),face_to_cell), face_to_cell,l)

q_Γ = GenericCellPoint(reindex(get_array(q),face_to_cell))
w_Γ = reindex(w,face_to_cell)
dΓ_ref = CellQuadrature(q_Γ,w_Γ)

ϕ_Γ_2 = ϕ∘refface_to_refcell_map
r1 = ϕ_Γ_1(q_Γ)
r2 = ϕ_Γ_2(q_Γ)
test_array(get_array(r1),get_array(r2))
test_array(get_cell_id(r1),get_cell_id(r2))

p = refface_to_refcell_map(q_Γ)
x = ϕ(p)
test_array(get_array(x),get_array(r2))
test_array(get_cell_id(x),get_cell_id(r2))

ϕ_Γ = FaceMap(ϕ_Γ_1,ϕ,refface_to_refcell_map)
test_cell_map(ϕ_Γ)

r1 = ϕ_Γ(q_Γ)
test_array(get_array(r1),get_array(r2))
test_array(get_cell_id(r1),get_cell_id(r2))

h_ref = convert_to_cell_field(x->x[1],length(face_to_cell))
h = h_ref∘inverse_map(ϕ_Γ)
@test (h∘ϕ_Γ)(q_Γ) == h_ref(q_Γ)

r1 = (dv∘ϕ_Γ)(q_Γ)
r2 = (s∘refface_to_refcell_map)(q_Γ)
test_array(get_array(r1),get_array(r2))

dΓ = ϕ_Γ(dΓ_ref)

r = ∫( 3*∇(dv)⋅∇(du)*f + h*dv*du*2 )*dΓ
test_array(r,collect(r))

r = ∫( 1 )*dΓ
test_array(r,collect(r))

r = ∫( x->1 )*dΓ
test_array(r,collect(r))

# Face map + re-indexed

ids = [2,1]
ϕ1_Γ = ReindexedCellMap(ϕ_Γ,ids)

q1_Γ = GenericCellPoint(reindex(get_array(q_Γ),ids))
w1_Γ = reindex(w_Γ,ids)
dΓ1_ref = CellQuadrature(q1_Γ,w1_Γ)

r1 = (dv∘ϕ1_Γ)(q1_Γ)
r2 = reindex( (dv∘ϕ_Γ)(q_Γ), ids)
test_array(r1,r2)

dΓ1 = ϕ1_Γ(dΓ1_ref)

r = ∫( 1 )*dΓ1
test_array(r,collect(r))

j_ref = convert_to_cell_field(x->x[1],length(ids))
j = j_ref∘inverse_map(ϕ1_Γ)
@test length(ϕ1_Γ) == 2
@test num_cell_ids(ϕ_Γ) == 4
@test num_cell_ids(ϕ1_Γ) == 4
@test length(j_ref) == 2
@test length(j) == 4

r = ∫( 3*∇(dv)⋅∇(du)*h + dv*du*j )*dΓ1
test_array(r,collect(r))

## Skeleton related

ϕ_Λ = SkeletonFaceMap(ϕ_Γ,ϕ_Γ)
test_cell_map(ϕ_Γ)
q_Λ = q_Γ

r = (dv.⁺∘ϕ_Λ)(q_Λ)
@test isa(r,VectorOfBlockArrayCoo)
test_array(r,collect(r))

r = (du.⁻∘ϕ_Λ)(q_Λ)
@test isa(r,VectorOfBlockArrayCoo)
test_array(r,collect(r))

r1 = (∇(dv).⁺∘ϕ_Λ)(q_Λ)
r2 = (∇(dv.⁺)∘ϕ_Λ)(q_Λ)
test_array(r1,r2)

r1 = (jump(∇(du))∘ϕ_Λ)(q_Λ)
r2 = (∇(jump(du))∘ϕ_Λ)(q_Λ)
test_array(r1,r2)

dΛ = ϕ_Λ(dΓ_ref)

r = ∫( h*mean(∇(dv))⋅jump(∇(du)) + du.⁺*jump(dv)  )*dΛ
@test isa(r,VectorOfBlockArrayCoo)
test_array(r,collect(r))

r = ∫( 1 )*dΛ
test_array(r,collect(r))

r = ∫( f )*dΛ
test_array(r,collect(r))

# ReindexedCellMap

ids1 = [4,2,3]
ϕ1 = ReindexedCellMap(ϕ,ids1)

q1 = GenericCellPoint(reindex(get_array(q),ids1))
w1 = reindex(w,ids1)
dΩ1_ref = CellQuadrature(q1,w1)

r1 = (dv∘ϕ1)(q1)
r2 = reindex( (dv∘ϕ)(q), ids1)
test_array(r1,r2)

dΩ1 = ϕ1(dΩ1_ref)

r = ∫( 1 )*dΩ1
test_array(r,collect(r))

h_ref = convert_to_cell_field(x->x[1],length(ids1))
h = h_ref∘inverse_map(ϕ1)
@test length(ϕ1) == 3
@test num_cell_ids(ϕ1) == 4
@test length(h_ref) == 3
@test length(h) == 4

r = ∫( 3*∇(dv)⋅∇(du)*h + dv*du*2 )*dΩ1
test_array(r,collect(r))

# lazy_append
ids2 = [3,1,4,2,3]
ϕ2 = ReindexedCellMap(ϕ,ids2)
ϕ12 = lazy_append(ϕ1,ϕ2)

q2 = GenericCellPoint(reindex(get_array(q),ids2))
w2 = reindex(w,ids2)
dΓ2_ref = CellQuadrature(q2,w2)

dΓ12_ref = lazy_append(dΩ1_ref,dΓ2_ref)
dΓ12 = ϕ12(dΓ12_ref)

q12 = get_coordinates(dΓ12_ref)
r = (dv∘ϕ12)(q12)
@test isa(r,AppendedArray)
test_array(r,collect(r))

r = ∫(1)*dΓ12
@test isa(r,AppendedArray)
test_array(r,collect(r))

r = ∫(x->1)*dΓ12
@test isa(r,AppendedArray)
test_array(r,collect(r))

r = ∫( 3*∇(dv)⋅∇(du)*f + 2*dv*du )*dΓ12
@test isa(r,AppendedArray)
test_array(r,collect(r))

n_ref(x) = x[1]
n1_ref = convert_to_cell_field(n_ref,length(ϕ1))
n2_ref = convert_to_cell_field(n_ref,length(ϕ2))
n12_ref = lazy_append(n1_ref,n2_ref)
n12 = n12_ref∘inverse_map(ϕ12)
r = (n12∘ϕ12)(q12)
@test isa(r,AppendedArray)
test_array(r,collect(r))

r = ∫( 3*∇(dv)⋅∇(du)*n12 + 2*dv*du )*dΓ12
@test isa(r,AppendedArray)
test_array(r,collect(r))

end # module
