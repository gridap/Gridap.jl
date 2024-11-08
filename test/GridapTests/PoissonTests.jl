module PoissonTests

using Test
using Gridap
using Gridap.TensorValues
import Gridap: ∇
#using LinearAlgebra

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

#u_scal(x) = x[1]^2 + x[2]
#∇u_scal(x) = VectorValue( 2*x[1], one(x[2]) )
#Δu_scal(x) = 2
#f_scal(x) = - Δu_scal(x)
#∇(::typeof(u_scal)) = ∇u_scal

scalar_data = Dict{Symbol,Any}()
scalar_data[:valuetype] = Float64
scalar_data[:u] = u_scal
scalar_data[:f] = f_scal

# Using automatic differentiation
u_vec(x) = VectorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2 )
f_vec(x) = - Δ(u_vec)(x)

#u_vec(x) = VectorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2 )
#∇u_vec(x) = TensorValue( 2*x[1], one(x[2]), 4*one(x[1]), - 2*x[2] )
#Δu_vec(x) = VectorValue( 2, -2 )
#f_vec(x) = - Δu_vec(x)
#∇(::typeof(u_vec)) = ∇u_vec

vector_data = Dict{Symbol,Any}()
vector_data[:valuetype] = VectorValue{2,Float64}
vector_data[:u] = u_vec
vector_data[:f] = f_vec

u_ten(x) = TensorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2, 2x[2]^2 - 3x[1], -.5x[1]^2 + x[2] )
f_ten(x) = - Δ(u_ten)(x)
tensor_data = Dict{Symbol,Any}()
tensor_data[:valuetype] = TensorValue{2,2,Float64}
tensor_data[:u] = u_ten
tensor_data[:f] = f_ten


u_sten(x) = SymTensorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2, 2x[2]^2 - 3x[1])
f_sten(x) = - Δ(u_sten)(x)
stensor_data = Dict{Symbol,Any}()
stensor_data[:valuetype] = SymTensorValue{2,Float64}
stensor_data[:u] = u_sten
stensor_data[:f] = f_sten


u_qten(x) = SymTracelessTensorValue( x[1]^2 + x[2], 4*x[1] - x[2]^2)
f_qten(x) = - Δ(u_qten)(x)
qtensor_data = Dict{Symbol,Any}()
qtensor_data[:valuetype] = SymTracelessTensorValue{2,Float64}
qtensor_data[:u] = u_qten
qtensor_data[:f] = f_qten

for data in [scalar_data, vector_data, tensor_data, stensor_data, qtensor_data]

  T = data[:valuetype]
  u = data[:u]
  f = data[:f]

  for domain_style in (ReferenceDomain(),PhysicalDomain())

    cell_fe = FiniteElements(domain_style,model,lagrangian,T,order)
    V = TestFESpace(Ω,cell_fe,dirichlet_tags="dirichlet",labels=labels)
    U = TrialFESpace(V,u)

    uh = interpolate(u, U)

    a(u,v) =
      ∫( ∇(v)⊙∇(u) )*dΩ +
      ∫( (γ/h)*v⊙u  - v⊙(nd⋅∇(u)) - (nd⋅∇(v))⊙u )*dΓd

    l(v) =
      ∫( v⊙f )*dΩ +
      ∫( v⊙(nn⋅∇(uh)) )*dΓn +
      ∫( (γ/h)*v⊙uh - (nd⋅∇(v))⊙uh )*dΓd

    op = AffineFEOperator(a,l,U,V)
    uh = solve(op)

    e = u - uh

    l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
    h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

    el2 = l2(e)
    eh1 = h1(e)
    ul2 = l2(uh)
    uh1 = h1(uh)

    @test el2/ul2 < 1.e-8
    @test eh1/uh1 < 1.e-7

  end
end

end # module
