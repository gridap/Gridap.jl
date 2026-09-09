module ArnoldWintherConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.Arrays, Gridap.TensorValues

# Hellinger--Reissner mixed elasticity. With the identity compliance, find
# (σ, u) ∈ Σₕ × Vₕ with
#
#   ∫ σ⊙τ + b(τ,u) = 0,   b(σ,v) = -∫ f⋅v,     b(τ,v) = Σ_K ∫_K div τ ⋅ v,
#
# and Vₕ the discontinuous vector P₁ space. The displacement boundary condition
# u = 0 is natural here, and the manufactured u vanishes on ∂Ω.
#
# `b` is integrated by parts element-wise, so that no derivative of τ is ever
# taken, and each side of a facet contributes with its own traction -- which
# keeps the form consistent for a nonconforming Σₕ. `solve_elasticity_div` below
# instead uses the `DoubleContraVariantPiolaMap` specialization of `DIV`; the two
# agree.

const π2 = pi^2
uex(x) = VectorValue(sin(pi * x[1]) * sin(pi * x[2]), sin(2*pi * x[1]) * sin(pi * x[2]))

function σex(x)
  u1x = pi * cos(pi * x[1]) * sin(pi * x[2])
  u1y = pi * sin(pi * x[1]) * cos(pi * x[2])
  u2x = 2*pi * cos(2*pi * x[1]) * sin(pi * x[2])
  u2y = pi * sin(2*pi * x[1]) * cos(pi * x[2])
  SymTensorValue{2,Float64}(u1x, (u1y + u2x) / 2, u2y)
end

function fex(x)
  u1 = sin(pi * x[1]) * sin(pi * x[2])
  u2 = sin(2*pi * x[1]) * sin(pi * x[2])
  d11u1 = -π2 * u1
  d22u1 = -π2 * u1
  d12u1 = π2 * cos(pi * x[1]) * cos(pi * x[2])
  d11u2 = -4*π2 * u2
  d22u2 = -π2 * u2
  d12u2 = 2*π2 * cos(2*pi * x[1]) * cos(pi * x[2])
  VectorValue(-(d11u1 + (d22u1 + d12u2) / 2), -((d12u1 + d11u2) / 2 + d22u2))
end

function solve_elasticity(nc, reffe_σ, degree; permute=false)
  model = simplexify(CartesianDiscreteModel((0,1,0,1), nc))
  permute && (model = permute_cells(model))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, degree)
  nΛ = get_normal_vector(Λ)
  Γ = BoundaryTriangulation(model)
  dΓ = Measure(Γ, degree)
  nΓ = get_normal_vector(Γ)

  Vσ = FESpace(model, reffe_σ)
  Vu = FESpace(model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1); conformity=:L2)
  X = MultiFieldFESpace([Vσ, Vu])
  Y = MultiFieldFESpace([Vσ, Vu])

  b(τ, v) = ∫( -(τ ⊙ ε(v)) )dΩ +
            ∫( ((τ.⁺ ⋅ nΛ.⁺) ⋅ v.⁺) + ((τ.⁻ ⋅ nΛ.⁻) ⋅ v.⁻) )dΛ +
            ∫( ((τ ⋅ nΓ) ⋅ v) )dΓ
  A((σ, u), (τ, v)) = ∫( σ ⊙ τ )dΩ + b(τ, u) + b(σ, v)
  L((τ, v)) = ∫( -(fex ⋅ v) )dΩ

  σh, uh = solve(AffineFEOperator(A, L, X, Y))
  (sqrt(sum(∫((σh - σex) ⊙ (σh - σex))dΩ)),
   sqrt(sum(∫((uh - uex) ⋅ (uh - uex))dΩ)))
end

# The same system with b(τ,v) written directly as ∫(DIV(τ)⋅v)dω against a
# reference-domain measure.
function solve_elasticity_div(nc, reffe_σ, degree)
  model = simplexify(CartesianDiscreteModel((0,1,0,1), nc))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  dω = Measure(Ω, degree, integration_domain_style=ReferenceDomain())

  Vσ = FESpace(model, reffe_σ)
  Vu = FESpace(model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1); conformity=:L2)
  X = MultiFieldFESpace([Vσ, Vu])
  Y = MultiFieldFESpace([Vσ, Vu])

  A((σ, u), (τ, v)) = ∫( σ ⊙ τ )dΩ + ∫( (DIV(τ) ⋅ u) + (DIV(σ) ⋅ v) )dω
  L((τ, v)) = ∫( -(fex ⋅ v) )dΩ

  σh, uh = solve(AffineFEOperator(A, L, X, Y))
  (sqrt(sum(∫((σh - σex) ⊙ (σh - σex))dΩ)),
   sqrt(sum(∫((uh - uex) ⋅ (uh - uex))dΩ)))
end

# The same mesh with the vertices of each cell listed in a permuted order: the
# discrete space must not depend on it.
function permute_cells(model, perms=([1,2,3],[2,3,1],[3,1,2],[2,1,3],[1,3,2],[3,2,1]))
  grid = get_grid(model)
  cell_nodes = get_cell_node_ids(grid)
  permuted = [collect(cell_nodes[i])[perms[mod1(i, length(perms))]]
              for i in 1:length(cell_nodes)]
  pgrid = UnstructuredGrid(get_node_coordinates(grid), Table(permuted),
                           get_reffes(grid), get_cell_type(grid))
  pmodel = UnstructuredDiscreteModel(pgrid)
  topo = get_grid_topology(pmodel)
  labels = get_face_labeling(pmodel)
  for d in 0:1
    get_face_entity(labels, d) .= ifelse.(Geometry.get_isboundary_face(topo, d), 2, 1)
  end
  get_face_entity(labels, 2) .= 1
  add_tag!(labels, "interior", [1])
  add_tag!(labels, "boundary", [2])
  pmodel
end

# The approximation power of the space alone, by interpolation
function interp_test(ncs, reffe)
  es = Float64[]
  for nc in ncs
    model = simplexify(CartesianDiscreteModel((0,1,0,1), nc))
    dΩ = Measure(Triangulation(model), 12)
    σh = interpolate(σex, FESpace(model, reffe))
    push!(es, sqrt(sum(∫((σh - σex) ⊙ (σh - σex))dΩ)))
  end
  es
end

function convg_test(ncs, reffe, degree)
  eσ = Float64[]
  eu = Float64[]
  hs = Float64[]
  for nc in ncs
    l2σ, l2u = solve_elasticity(nc, reffe, degree)
    println((l2σ, l2u))
    push!(eσ, l2σ)
    push!(eu, l2u)
    h = 1/nc[1]
    push!(hs, h)
  end
  println(eσ)
  println(eu)
  println(hs)
  eσ, eu, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

ncs = [(4,4),(8,8),(16,16)]

# AWnc: the space approximates σ at O(h²), while the nonconforming method
# delivers O(h) in the stress and O(h²) in the displacement -- the low end of the
# published range, not a defect.
reffe_nc = ReferenceFE(TRI, aw_nc, Float64)
println("--- AWnc ---")
println("Interpolation of σ: $(interp_test(ncs, reffe_nc))")
eσ, eu, hs = convg_test(ncs, reffe_nc, 8)
println("Slope L2-norm σ: $(slope(hs,eσ))")
println("Slope L2-norm u: $(slope(hs,eu))")
println("Sorted vs permuted: $(solve_elasticity((8,8), reffe_nc, 8)) vs $(solve_elasticity((8,8), reffe_nc, 8; permute=true))")
println("By parts vs DIV: $(solve_elasticity((4,4), reffe_nc, 8)) vs $(solve_elasticity_div((4,4), reffe_nc, 8))")

# AWc: O(h³) in the stress and O(h²) in the displacement
reffe_c = ReferenceFE(TRI, aw_c, Float64)
println("--- AWc ---")
eσ, eu, hs = convg_test(ncs, reffe_c, 10)
println("Slope L2-norm σ: $(slope(hs,eσ))")
println("Slope L2-norm u: $(slope(hs,eu))")
println("Sorted vs permuted: $(solve_elasticity((8,8), reffe_c, 10)) vs $(solve_elasticity((8,8), reffe_c, 10; permute=true))")
println("By parts vs DIV: $(solve_elasticity((4,4), reffe_c, 10)) vs $(solve_elasticity_div((4,4), reffe_c, 10))")

end # module
