module AssemblersTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.Algebra
using SparseArrays
using Gridap.FESpaces
using Gridap.CellData

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

V = TestFESpace(
  model,
  ReferenceFE(lagrangian,Float64,1),
  dirichlet_tags=[1,2,3,4,6,5])

u(x) = x[1]+x[2]

U = TrialFESpace(u,V)

dv = get_fe_basis(V)
du = get_trial_fe_basis(U)

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,degree)

a(u,v) = ∫(∇(u)⋅∇(v))*dΩ + ∫(u*v)*dΓ
ℓ(v) = ∫(v)*dΩ

mat_contribs = a(du,dv)
vec_contribs = ℓ(dv)

assem = SparseMatrixAssembler(U,V)

@test isa(U,TrialFESpace)
#@test_throws AssertionError assem = SparseMatrixAssembler(V,U)

data = collect_cell_matrix(U,V,mat_contribs)
A = assemble_matrix(assem,data)
@test size(A) == (num_free_dofs(V), num_free_dofs(U))

data = collect_cell_vector(V,vec_contribs)
b = assemble_vector(assem,data)
x = A\b
uh = FEFunction(U,x)

data = collect_cell_matrix_and_vector(U,V,mat_contribs,vec_contribs)
A,b = assemble_matrix_and_vector(assem,data)
x = A\b
uh = FEFunction(U,x)

uhd = zero(U)
data = collect_cell_matrix_and_vector(U,V,mat_contribs,vec_contribs,uhd)
A,b = assemble_matrix_and_vector(assem,data)
x = A\b
uh = FEFunction(U,x)

A1 = assemble_matrix(a,U,V)
b1 = assemble_vector(ℓ,V)
A2,b2 = assemble_matrix_and_vector(a,ℓ,U,V)

A12 = copy(A1); A12[1,1]=rand()
b12 = copy(b1); b12[1]=rand()
A22 = copy(A2); A22[1,1]=rand()
b22 = copy(b2); b22[1]=rand()

tol = 1.e-14
assemble_matrix!(a,A12,U,V)
assemble_vector!(ℓ,b12,V)
assemble_matrix_and_vector!(a,ℓ,A22,b22,U,V)
@test norm(A12-A1) < tol
@test norm(b12-b1) < tol
@test norm(A22-A2) < tol
@test norm(b22-b2) < tol

A1 = assemble_matrix(a(du,dv),U,V)
b1 = assemble_vector(ℓ(dv),V)
A2,b2 = assemble_matrix_and_vector(a(du,dv),ℓ(dv),U,V)

A12 = copy(A1); A12[1,1]=rand()
b12 = copy(b1); b12[1]=rand()
A22 = copy(A2); A22[1,1]=rand()
b22 = copy(b2); b22[1]=rand()

assemble_matrix!(a(du,dv),A12,U,V)
assemble_vector!(ℓ(dv),b12,V)
assemble_matrix_and_vector!(a(du,dv),ℓ(dv),A22,b22,U,V)
@test norm(A12-A1) < tol
@test norm(b12-b1) < tol
@test norm(A22-A2) < tol
@test norm(b22-b2) < tol

# Rectangular matrix
a2(u,v) = ∫(u⋅∇(v))*dΩ
V2 = TestFESpace(
  model,
  ReferenceFE(lagrangian,VectorValue{2,Float64},1))
U2 = TrialFESpace(V2)

du2 = get_trial_fe_basis(U2)
mat_contribs = a2(du2,dv)
vec_contribs = ℓ(dv)

assem = SparseMatrixAssembler(U2,V)
data = collect_cell_matrix(U2,V,mat_contribs)
A2 = assemble_matrix(assem,data)
@test size(A2) == (num_free_dofs(V), num_free_dofs(U2))
A2,b2 = assemble_matrix_and_vector(a2(du2,dv),ℓ(dv),U2,V)
@test size(A2) == (num_free_dofs(V), num_free_dofs(U2))
A2,b2 = assemble_matrix_and_vector(a2(du2,dv),0,U2,V)
@test size(A2) == (num_free_dofs(V), num_free_dofs(U2))

V = TestFESpace(
  model,
  ReferenceFE(lagrangian,Float64,1),
  vector_type=Vector{ComplexF64})
U = V

assem = SparseMatrixAssembler(U,V)
Ta = get_matrix_type(assem)
Tb = get_vector_type(assem)
@test eltype(Ta) == ComplexF64
@test eltype(Tb) == ComplexF64

# Now with an homogeneous linear form

a(u,v) = ∫(∇(u)⋅∇(v))*dΩ + ∫(u*v)*dΓ
ℓ(v) = 0

mat_contribs = a(du,dv)
vec_contribs = ℓ(dv)

data = collect_cell_matrix(U,V,mat_contribs)
A = assemble_matrix(assem,data)
@test size(A) == (num_free_dofs(V), num_free_dofs(U))

data = collect_cell_vector(V,vec_contribs)
b = assemble_vector(assem,data)
x = A\b
@test x ≈ b

data = collect_cell_matrix_and_vector(U,V,mat_contribs,vec_contribs)
A,b = assemble_matrix_and_vector(assem,data)
x = A\b
@test x ≈ b

uhd = zero(U)
data = collect_cell_matrix_and_vector(U,V,mat_contribs,vec_contribs,uhd)
A,b = assemble_matrix_and_vector(assem,data)
x = A\b
@test x ≈ b


mutable struct MyArrayCounter
  count
  MyArrayCounter()=new(0)
end
function Algebra.LoopStyle(::Type{<:MyArrayCounter})
  Loop()
end

function Algebra.add_entry!(c::Function,a::MyArrayCounter,v,i,j)
  a.count=a.count+1
end
function Algebra.add_entry!(c::Function,a::MyArrayCounter,v,i)
  a.count=a.count+1
end
function Algebra.add_entries!(c::Function,a::MyArrayCounter,v,i,j)
  a.count=a.count+length(i)
end
function Algebra.add_entries!(c::Function,a::MyArrayCounter,v,i)
  a.count=a.count+length(i)
end

mac=MyArrayCounter()
vec_contribs = ∫(1*dv)*dΩ
data = collect_cell_vector(V,vec_contribs)
symbolic_loop_vector!(mac,assem,data)
@test mac.count == num_cells(Ω)*4




end # module
