module MomentBasedReferenceFEsTests

using Test
using Gridap
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs

# Tests for differential operators
# We assume that taking integral just work, and only test taking point values
T = Float64

p = TET
D = 3
V = VectorValue{D,T}

vertices_range = get_dimrange(p,0)
poly_vertices = get_vertex_coordinates(p)

function comps_mom(Dφ,μ,ds)
  Broadcasting(Operation(⊙))(Dφ,μ)
end


φ = MonomialBasis(Val(D),V,2)
μ = MonomialBasis(Val(0),V,0)
V_moments = [ (vertices_range, comps_mom, μ), ]
dofs_V = MomentBasedDofBasis(p,φ,V_moments)
Vφv = reinterpret(T, evaluate(φ, poly_vertices))
@test Vφv == evaluate(dofs_V, φ)

Gφ = Broadcasting(∇)(φ)
G = gradient_type(V, zero(Point{D,T}))
μG = MonomialBasis(Val(0),G,0)
G_moments = [ (vertices_range, comps_mom, μG), ]
dofs_G = MomentBasedDofBasis(p,φ,G_moments,∇)
Gφv = reinterpret(T, evaluate(Gφ, poly_vertices))
@test Gφv == evaluate(dofs_G, φ)

Hφ = Broadcasting(∇)(Gφ)
H = gradient_type(G, zero(Point{D,T}))
μH = MonomialBasis(Val(0),H,0)
H_moments = [ (vertices_range, comps_mom, μH), ]
dofs_H = MomentBasedDofBasis(p,φ,H_moments,∇∇)
Hφv = reinterpret(T, evaluate(Hφ, poly_vertices))
@test Hφv == evaluate(dofs_H, φ)

Dφ = Broadcasting(divergence)(φ)
μD = MonomialBasis(Val(0),T,0) # div on 3D-vector results in scalar
D_moments = [ (vertices_range, comps_mom, μD), ]
dofs_D = MomentBasedDofBasis(p,φ,D_moments,divergence)
Dφv = reinterpret(T, evaluate(Dφ, poly_vertices))
@test Dφv == evaluate(dofs_D, φ)

Cφ = Broadcasting(curl)(φ)
μC = μ # curl on 3D-vector results in 3D-vectors
C_moments = [ (vertices_range, comps_mom, μC), ]
dofs_C = MomentBasedDofBasis(p,φ,C_moments,curl)
Cφv = reinterpret(T, evaluate(Cφ, poly_vertices))
@test Cφv == evaluate(dofs_C, φ)


p = TRI
D = 2
P = Point{D,T}
for V in (SymTracelessTensorValue{3,T}, SkewSymTensorValue{3,T})
  φ = MonomialBasis(Val(D),V,2)
  Gφ = Broadcasting(∇)(φ)
  G = gradient_type(V, zero(P))
  grad_V_basis = [ p⊗v for p in component_basis(P) for v in component_basis(V)]
  grad_V_dual_basis = representatives_of_basis_dual(grad_V_basis)

  # Create a polynomial basis having constant basis polynomial equal to grad_V_basis
  e = MonomialBasis(Val(0),G,0)    # cartesian vector basis
  e_basis = evaluate(e, Point())   # actually same as component_basis(G) (or permutation of)
  change = [ gd⊙eᵢ for eᵢ in e_basis, gd in grad_V_dual_basis ]
  μ = linear_combination(change, e)
  @test grad_V_dual_basis == evaluate(μ, Point())

  moments = [ (get_dimrange(p,0), comps_mom, μ), ]

  dofs_G = MomentBasedDofBasis(p,φ,moments,∇)
  dofs_G_φ = evaluate(dofs_G, φ)
  poly_vertices = get_vertex_coordinates(p)
  Gφv = evaluate(Gφ, poly_vertices)
  # compute the gradients with respect to independat components of Gφ by hand
  Gφiv = T[ Gφv[k,j]⊙G_dual_i for G_dual_i in grad_V_dual_basis, k in 1:size(Gφv,1), j in 1:size(Gφv,2)]
  Gφiv = reshape(Gφiv, size(dofs_G_φ))

  @test Gφiv == dofs_G_φ
end

end
