
using Test
using Gridap
using Gridap.ReferenceFEs, Gridap.Geometry
using Gridap.Helpers
using Gridap.FESpaces
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Fields

using LinearAlgebra


function nonpermutedmodel()
  conn   = Table([[1 2 3],[4,3,2]])
  coords = [VectorValue(0.0, 0.0),VectorValue(1.0, 0.0),VectorValue(0.0, 1.0),VectorValue(1.0, 1.0)]
  polys  = [TRI]
  reffes = [TRI3]
  ctypes = Int8[1,1]

  topo = Geometry.UnstructuredGridTopology(coords,conn,ctypes,polys,Geometry.NonOriented())
  grid = Geometry.UnstructuredGrid(coords,conn,reffes,ctypes,Geometry.NonOriented())
  labels = Geometry.FaceLabeling(topo)
  model = UnstructuredDiscreteModel(grid,topo,labels)
end

using Gridap.CellData: SkeletonPair
function LinearAlgebra.dot(a::SkeletonPair,b::SkeletonPair)
  SkeletonPair(dot(a.plus,b.plus),dot(a.minus,b.minus))
end

function l2_error(uh,u,dΩ)
  eh = uh - u
  println(sum(∫(eh⊙eh)dΩ))
  println(collect(get_array(∫(eh⊙eh)dΩ)))
end

function airy(f,x)
  H = ∇(∇(f))(x)
  SymTensorValue(H[2,2],-H[2,1],H[1,1])
end
airy(f) = x -> airy(f,x)

ϕ(x) = x[1]^3 - 3.0*x[2]^2 + 2.0*x[1]*x[2]
u(x) = airy(ϕ)(x)

p = TRI
D = num_dims(p)
order = 0

reffe = ReferenceFEs.HellanHerrmannJhonsonRefFE(Float64,p,order)
dofs_ref = get_dof_basis(reffe)
basis_ref = get_shapefuns(reffe) 
prebasis  = get_prebasis(reffe)
get_face_own_dofs(reffe)
get_face_own_dofs_permutations(reffe)

M = evaluate(dofs_ref, basis_ref)
@test M ≈ I(num_dofs(reffe))

model = simplexify(CartesianDiscreteModel((0,1,0,1), (1,1));positive=true)
#model = nonpermutedmodel()
topo = get_grid_topology(model)
cell_to_nodes = Geometry.get_faces(topo,2,0)
cell_to_faces = Geometry.get_faces(topo,2,1)
face_to_nodes = Geometry.get_faces(topo,1,0)
face_to_cells = Geometry.get_faces(topo,1,2)
perms = get_cell_permutations(topo)
cmaps = get_cell_map(get_grid(model))
Jt = lazy_map(Broadcasting(∇),cmaps)

pts = [VectorValue(0.,0.),VectorValue(1.0,0.),VectorValue(0.,1.)]

Ω = Triangulation(model)
Γ = Boundary(model)
Γ.glue.face_to_bgface
Λ = Skeleton(model)
Λ.plus.glue.face_to_lcell
Λ.minus.glue.face_to_lcell

dΩ = Measure(Ω,4*(order+1))
dΓ = Measure(Γ,4*(order+1))
dΛ = Measure(Λ,4*(order+1))

V = FESpace(model,reffe)
get_cell_dof_ids(V)

uh = interpolate(u,V)

x = zeros(num_free_dofs(V))
#x[[1,2,4,5]] .= 1.0
x[3] = 1.0
uh = FEFunction(V,x)
get_cell_dof_values(uh)

uh_bis = interpolate(uh,V)
get_cell_dof_values(uh_bis)

cf = uh
get_free_dof_values(cf)

print_op_tree(cf.cell_field.cell_field)

#μ_c = MonomialBasis(Val(2),SymTensorValue{2,Float64},order-1,Polynomials._p_filter)
#μ_c = get_shapefuns(LagrangianRefFE(SymTensorValue{2,Float64}, TRI, order-1))
μ_c = map(constant_field,[SymTensorValue(0.,1.,0.),SymTensorValue(-2.,1.,0.),SymTensorValue(0.,-1.,2.)])
μ_Ω = GenericCellField([μ_c,μ_c],Ω,ReferenceDomain())
collect(get_array(∫((cf⊙transpose(μ_Ω)))dΩ))

# μ_f = MonomialBasis(Val(1),Float64,order,Polynomials._p_filter)
μ_f = get_shapefuns(LagrangianRefFE(Float64, SEGMENT, order))
μ_Γ = GenericCellField([μ_f,μ_f,μ_f,μ_f],Γ,ReferenceDomain())
μ_Λ = GenericCellField([μ_f],Λ,ReferenceDomain())
nΓ = get_normal_vector(Γ)
nΛ = get_normal_vector(Λ)

collect(get_array(∫(nΓ⋅(cf⋅nΓ)*μ_Γ)dΓ))
collect(get_array(∫((nΛ⋅(cf⋅nΛ)).plus - (nΛ⋅(cf⋅nΛ)).minus)dΛ))
collect(get_array(∫((nΛ⋅(cf⋅nΛ)).plus)dΛ))
collect(get_array(∫((nΛ⋅(cf⋅nΛ)).minus)dΛ))

pts_Ω = get_cell_points(dΩ.quad)
pts_Γ = get_cell_points(dΓ.quad)
pts_Λ = get_cell_points(dΛ.quad)

cf_Λ = cf#(cf⋅nΛ)
cf_Λ_plus = nΛ.plus⋅(cf.plus⋅nΛ.plus)
cf_Λ_minus = nΛ.minus⋅(cf.minus⋅nΛ.minus)

cf_plus = cf.plus
cf_minus = cf.minus

cf_Λ.plus(pts_Λ)
cf_Λ.minus(pts_Λ)

cf_Λ_plus(pts_Λ)
cf_Λ_minus(pts_Λ)

cf_plus(pts_Λ)[1]
cf_minus(pts_Λ)[1]

cf_plus_Λ = change_domain(cf_plus,Λ,ReferenceDomain())
print_op_tree(get_data(cf_plus_Λ))
cf_plus_Λ(pts_Λ)


pts = get_cell_points(dΩ.quad)
dofs_phys = get_fe_dof_basis(V)
basis_phys = get_fe_basis(V)
map(a -> round.(a),collect(dofs_phys(basis_phys)))[2]
y = collect(basis_phys(pts))

reffe_lag = LagrangianRefFE(SymTensorValue{2,Float64},TRI,1)
