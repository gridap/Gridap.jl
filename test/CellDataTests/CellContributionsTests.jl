module CellContributionsTests

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
dΩ_ref = CellQuadrature(q,w)
dΩ = LebesgueMeasure(ϕ(dΩ_ref))

# FaceMap
cell_map = ϕ
face_to_cell = [2,3,1]
cf = convert_to_cell_field(identity,length(face_to_cell))
refface_to_refcell_map = GenericCellMap(get_array(cf),face_to_cell,l)
ϕ_Γ_1 = GenericCellMap(reindex(get_array(ϕ),face_to_cell), face_to_cell,l)
ϕ_Γ = FaceMap(ϕ_Γ_1,ϕ,refface_to_refcell_map)

q_Γ = GenericCellPoint(reindex(get_array(q),face_to_cell))
w_Γ = reindex(w,face_to_cell)
dΓ_ref = CellQuadrature(q_Γ,w_Γ)
dΓ = LebesgueMeasure(ϕ_Γ(dΓ_ref))

# Functions
# TODO For the moment, we need this.
s_with_physical_grad = GenericCellField(attachmap(get_array(s),get_array(ϕ)),get_cell_axes(s),MetaSizeStyle(s))
dv = s_with_physical_grad∘inverse_map(ϕ)
du = trialize_cell_basis(dv)
u(x) = x[1]+x[2]
cell_to_u = [ u.(node_to_x[nodes]) for nodes in cell_to_node]
uh = lincomb(dv,cell_to_u)

a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ + ∫( u*v )*dΓ
b(v) = ∫( ∇(v)⋅∇(uh) )*dΩ + ∫( uh*v )*dΓ

mat_terms = a(du,dv)
vec_terms = b(dv)

a(u,v) = ∫( u*v )*dΓ
b(v) = ∫( uh*v )*dΓ

mat_terms = a(du,dv)
vec_terms = b(dv)

test_array(mat_terms,get_array(mat_terms))
test_array(vec_terms,get_array(vec_terms))

#matdata = collect_cell_matrix(mat_terms)
#vecdata = collect_cell_vector(vec_terms)
#matvecdata, matdata, vecdata = collect_cell_matrix_and_vector(mat_terms,vec_terms)
#matvecdata, matdata, vecdata = collect_cell_matrix_and_vector(mat_terms,vec_terms,uh)

end # module
