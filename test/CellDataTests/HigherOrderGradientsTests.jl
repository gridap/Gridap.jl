module HigherOrderGradientsTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using LinearAlgebra

# Polynomials interpolated exactly in the FE space, on Cartesian meshes with non-square cells
# (so the affine pushforward is not the identity), against their derivatives in closed form.
# A vector polynomial is a tuple of components, each a list of (exponents => coefficient).

function dmonomial(α::NTuple{D,Int}, dirs, x) where D
  β = collect(α)
  c = 1.0
  for d in dirs
    β[d] == 0 && return 0.0
    c *= β[d]
    β[d] -= 1
  end
  c * prod(x[d]^β[d] for d in 1:D)
end

evalpoly(P, x, dirs = ()) = sum(c * dmonomial(α, dirs, x) for (α, c) in P; init = 0.0)

# Entries of the exact N-th gradient in Gridap's storage order: derivative indices first
# (fastest), then the value components.
exact_entries(Ps, x, N, D) = [evalpoly(P, x, Tuple(I)) for P in Ps for I in CartesianIndices(ntuple(_ -> D, N))]

entries(v::Number) = [v]
entries(v::MultiValue) = collect(Tuple(v))

function max_error(cf, exact, dΩ)
  q = get_cell_points(dΩ)
  vals = collect(cf(q))
  xs = collect(q.cell_phys_point)
  maximum(maximum(maximum(abs, entries(v) .- exact(x)) for (v, x) in zip(vc, xc)) for (vc, xc) in zip(vals, xs))
end

# Lagrangian Q2, 2D: every order up to the total degree, and zero above it

model = CartesianDiscreteModel((0,2,0,1), (4,3))
Ω = Triangulation(model)
dΩ = Measure(Ω, 6)
P = [(2,2)=>1.0, (2,1)=>-2.0, (1,2)=>0.5, (2,0)=>1.5, (0,1)=>1.0]
V = FESpace(model, ReferenceFE(lagrangian, Float64, 2))
uh = interpolate(x -> evalpoly(P, x), V)
for N in 1:5
  @test max_error(gradient(uh, Val(N)), x -> exact_entries((P,), x, N, 2), dΩ) < 1e-9
end
q = get_cell_points(dΩ)
@test all(collect(gradient(uh, Val(1))(q)) .== collect(∇(uh)(q)))
@test all(collect(gradient(uh, Val(2))(q)) .== collect(∇∇(uh)(q)))
@test all(all(iszero, v) for v in collect(gradient(uh, Val(5))(q)))

# Lagrangian Q2, 3D

model3 = CartesianDiscreteModel((0,1,0,2,0,1), (2,3,2))
dΩ3 = Measure(Triangulation(model3), 4)
P3 = [(2,2,2)=>1.0, (2,1,0)=>-1.0, (0,1,2)=>0.5]
V3 = FESpace(model3, ReferenceFE(lagrangian, Float64, 2))
uh3 = interpolate(x -> evalpoly(P3, x), V3)
for N in 1:3
  @test max_error(gradient(uh3, Val(N)), x -> exact_entries((P3,), x, N, 3), dΩ3) < 1e-10
end

# Piola-mapped bases, 2D: a Q1 vector polynomial lies in RT1 and ND1

Pv = ([(1,1)=>1.0, (1,0)=>0.5, (0,1)=>-1.0], [(1,1)=>-0.7, (0,1)=>1.2, (0,0)=>1.0])
for (reffe, conformity) in ((raviart_thomas, :Hdiv), (nedelec, :Hcurl))
  Vv = FESpace(model, ReferenceFE(reffe, Float64, 1); conformity)
  uv = interpolate(x -> VectorValue(evalpoly(Pv[1], x), evalpoly(Pv[2], x)), Vv)
  for N in 1:3
    @test max_error(gradient(uv, Val(N)), x -> exact_entries(Pv, x, N, 2), dΩ) < 1e-9
  end
end

# N-th gradients of FE bases assemble, also through jumps on the skeleton

u = get_trial_fe_basis(V)
v = get_fe_basis(V)
A = assemble_matrix(∫(gradient(u, Val(3)) ⊙ gradient(v, Val(3)))dΩ, V, V)
@test A ≈ A'
x = get_free_dof_values(uh)
@test x' * A * x ≈ sum(∫(gradient(uh, Val(3)) ⊙ gradient(uh, Val(3)))dΩ)

Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ, 6)
B = assemble_matrix(∫(jump(gradient(u, Val(3))) ⊙ jump(gradient(v, Val(3))))dΛ, V, V)
@test norm(B) > 0
@test B ≈ B'
@test abs(x' * B * x) < 1e-12 * norm(B) * norm(x)^2   # a global polynomial has no jumps

end # module
