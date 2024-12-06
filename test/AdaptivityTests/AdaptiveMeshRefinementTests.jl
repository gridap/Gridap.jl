using Test
using Gridap, Gridap.Geometry, Gridap.Adaptivity
using DataStructures

# Marking tests

function test_dorfler_marking()
  for strategy in (:sort,:binsort,:quickmark)
    for n in (1000,10000)
      for θ in (0.3,0.5)
        η = abs.(randn(n))
        m = DorflerMarking(θ;strategy)
        I = Adaptivity.mark(m,η)
        @test sum(η[I]) > θ * sum(η)
        println("Strategy: $strategy, θ: $θ, n: $n")
        println("  > N marked = $(length(I)), val marked = $(sum(η[I]) / sum(η))")
      end
    end
  end
end

@testset "Dorfler marking" test_dorfler_marking()

# AMR tests

function LShapedModel(n)
  model = CartesianDiscreteModel((0,1,0,1),(n,n))
  cell_coords = map(mean,get_cell_coordinates(model))
  l_shape_filter(x) = (x[1] < 0.5) || (x[2] < 0.5)
  mask = map(l_shape_filter,cell_coords)
  return simplexify(DiscreteModelPortion(model,mask))
end

function amr_step(model)
  order = 1
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe)
  
  Ω = Triangulation(model)
  Γ = Boundary(model)
  Λ = Skeleton(model)
  
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)
  dΛ = Measure(Λ,2*order)
  
  hK = CellField(sqrt.(collect(get_array(∫(1)dΩ))),Ω)
  
  f(x) = 1.0 / ((x[1]-0.5)^2 + (x[2]-0.5)^2)^(1/2)
  a(u,v) = ∫(∇(u)⋅∇(v))dΩ
  l(v)   = ∫(f*v)dΩ
  ηh(u)  = ∫(hK*f)dΩ + ∫(hK*∇(u)⋅∇(u))dΓ + ∫(hK*jump(∇(u))⋅jump(∇(u)))dΛ
  
  op = AffineFEOperator(a,l,V,V)
  uh = solve(op)
  η = estimate(ηh,uh)
  
  m = DorflerMarking(0.5)
  I = Adaptivity.mark(m,η)
  
  method = Adaptivity.NVBRefinement(model)
  fmodel = Adaptivity.get_model(refine(method,model;cells_to_refine=I))

  return fmodel, uh, η, I
end

model = LShapedModel(10)

for i in 1:10
  fmodel, uh, η, I = amr_step(model)
  is_refined = map(i -> ifelse(i ∈ I, 1, -1), 1:num_cells(model))
  Ω = Triangulation(model)
  writevtk(
    Ω,"tmp/model_$(i-1)",append=false,
    cellfields = [
      "uh" => uh,
      "η" => CellField(η,Ω),
      "is_refined" => CellField(is_refined,Ω)
    ],
  )
  model = fmodel
end
