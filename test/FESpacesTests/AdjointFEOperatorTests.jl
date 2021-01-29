module AdjointFEOperatorTests

using Test
using Gridap
using Gridap.FESpaces
using Gridap.Algebra
using SparseArrays: sparse

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)
order = 2

const h = (domain[2]-domain[1]) / partition[1]
const γ = 10

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[7,8])
add_tag_from_tags!(labels,"nitsche",6)

Ω = Triangulation(model)
Γn = BoundaryTriangulation(model,labels,tags="neumann")
Γd = BoundaryTriangulation(model,labels,tags="nitsche")

degree = order
dΩ = Measure(Ω,degree)
dΓn = Measure(Γn,degree)
dΓd = Measure(Γd,degree)

nn = get_normal_vector(Γn)
nd = get_normal_vector(Γd)

# Using automatic differentiation
u_scal(x) = x[1]^2 + x[2]
f_scal(x) = - Δ(u_scal)(x)

scalar_data = Dict{Symbol,Any}()
scalar_data[:valuetype] = Float64
scalar_data[:u] = u_scal
scalar_data[:f] = f_scal

# Using automatic differentiation
u_vec(x) = VectorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2 )
f_vec(x) = - Δ(u_vec)(x)

vector_data = Dict{Symbol,Any}()
vector_data[:valuetype] = VectorValue{2,Float64}
vector_data[:u] = u_vec
vector_data[:f] = f_vec

for data in [ vector_data, scalar_data ]

  T = data[:valuetype]
  u = data[:u]
  f = data[:f]

  for domain_style in (ReferenceDomain(),PhysicalDomain())

    cell_fe = FiniteElements(domain_style,model,lagrangian,T,order)
    V = TestFESpace(model,cell_fe,dirichlet_tags="dirichlet",labels=labels)
    U = TrialFESpace(V,u)

    uh = interpolate(u, U)

    a(u,v) =
      ∫( ∇(v)⊙∇(u) )*dΩ +
      ∫( (γ/h)*v⊙u  - v⊙(nd⋅∇(u)) - (nd⋅∇(v))⊙u )*dΓd

    l(v) =
      ∫( v⊙f )*dΩ +
      ∫( v⊙(nn⋅∇(uh)) )*dΓn +
      ∫( (γ/h)*v⊙uh - (nd⋅∇(v))⊙u )*dΓd

    ld(u) = ∫(u⋅u)dΩ

    op = AffineFEOperatorWithAdjoint(a,l,ld,U,V)
    op_forward = get_forward(op)
    op_adjoint = get_adjoint(op)

    uh,ψh = solve(op)
    uh_forward = solve(op_forward)
    ψh_adjoint = solve(op_adjoint)

    cellmat = get_matrix(op)
    cellmat_adj = get_adjoint_matrix(op)

    @test cellmat_adj == sparse(transpose(cellmat))

    l2(w) = √(∑(∫(w⊙w)dΩ))
    @time l2uh = l2(uh)
    @time l2ψh = l2(ψh)
    @time l2uh_forward = l2(uh_forward)
    @time l2ψh_adjoint = l2(ψh_adjoint)
    @test l2uh ≈ l2uh_forward
    @test l2ψh ≈ l2ψh_adjoint


    #res(u,v) = a(u,v) - l(v)
    #op_AD = FEOperator(res,U,V)
    #uh_AD = solve(op_AD)
    #adj_op = AdjointFEOperator(op,uh)
    #adj_op_AD = AdjointFEOperator(op_AD,uh_AD)
    #cellmat_adj_AD = get_matrix(adj_op_AD)
    #@test cellmat_adj ≈ cellmat_adj_AD
  end
end

end # module
