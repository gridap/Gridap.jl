module AssemblersTests

using Test

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
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
  ReferenceFE(:Lagrangian,Float64,1),
  dirichlet_tags=[1,2,3,4,6,5])

u(x) = x[1]+x[2]

U = TrialFESpace(u,V)

dv = get_cell_shapefuns(V)
du = get_cell_shapefuns_trial(U)

degree = 2
Ω = Triangulation(model)
dΩ = LebesgueMeasure(Ω,degree)

Γ = BoundaryTriangulation(model)
dΓ = LebesgueMeasure(Γ,degree)

a(u,v) = ∫(∇(u)⋅∇(v))*dΩ + ∫(u*v)*dΓ
ℓ(v) = ∫(v)*dΩ

mat_contribs = a(du,dv)
vec_contribs = ℓ(dv)

assem = SparseMatrixAssembler(U,V)

@test isa(U,TrialFESpace)
@test_throws AssertionError assem = SparseMatrixAssembler(V,U)

data = collect_cell_matrix(mat_contribs)
A = assemble_matrix(assem,data)
@test size(A) == (num_free_dofs(V), num_free_dofs(U))

data = collect_cell_vector(vec_contribs)
b = assemble_vector(assem,data)
x = A\b
uh = FEFunction(U,x)

data = collect_cell_matrix_and_vector(mat_contribs,vec_contribs)
A,b = assemble_matrix_and_vector(assem,data)
x = A\b
uh = FEFunction(U,x)

uhd = zero(U)
data = collect_cell_matrix_and_vector(mat_contribs,vec_contribs,uhd)
A,b = assemble_matrix_and_vector(assem,data)
x = A\b
uh = FEFunction(U,x)

A = assemble_matrix(a,U,V)
b = assemble_vector(ℓ,V)
A,b = assemble_matrix_and_vector(a,ℓ,U,V)

A = assemble_matrix(a(du,dv),U,V)
b = assemble_vector(ℓ(dv),V)
A,b = assemble_matrix_and_vector(a(du,dv),ℓ(dv),U,V)

V = TestFESpace(
  model,
  ReferenceFE(:Lagrangian,Float64,1),
  vector_type=Vector{ComplexF64})
U = V

assem = SparseMatrixAssembler(U,V)
Ta = get_matrix_type(assem)
Tb = get_vector_type(assem)
@test eltype(Ta) == ComplexF64
@test eltype(Tb) == ComplexF64

end # module
