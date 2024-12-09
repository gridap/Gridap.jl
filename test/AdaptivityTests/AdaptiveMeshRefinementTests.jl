module AdaptiveMeshRefinementTests

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

# AMR tests

function LShapedModel(n)
  model = CartesianDiscreteModel((0,1,0,1),(n,n))
  cell_coords = map(mean,get_cell_coordinates(model))
  l_shape_filter(x) = (x[1] < 0.5) || (x[2] < 0.5)
  mask = map(l_shape_filter,cell_coords)
  return simplexify(DiscreteModelPortion(model,mask))
end

l2_norm(he,xh,dΩ) = ∫(he*(xh*xh))*dΩ
l2_norm(xh,dΩ) = ∫(xh*xh)*dΩ

function amr_step(model,u_exact;order=1)
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe;dirichlet_tags=["boundary"])
  U = TrialFESpace(V,u_exact)
  
  Ω = Triangulation(model)
  Γ = Boundary(model)
  Λ = Skeleton(model)
  
  dΩ = Measure(Ω,4*order)
  dΓ = Measure(Γ,2*order)
  dΛ = Measure(Λ,2*order)
  
  hK = CellField(sqrt.(collect(get_array(∫(1)dΩ))),Ω)

  nΓ = get_normal_vector(Γ)
  nΛ = get_normal_vector(Λ)

  ∇u(x)  = ∇(u_exact)(x)
  f(x)   = -Δ(u_exact)(x)
  a(u,v) = ∫(∇(u)⋅∇(v))dΩ
  l(v)   = ∫(f*v)dΩ
  ηh(u)  = l2_norm(hK*(f + Δ(u)),dΩ) + l2_norm(hK*(∇(u) - ∇u)⋅nΓ,dΓ) + l2_norm(jump(hK*∇(u)⋅nΛ),dΛ)
  
  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)
  η = estimate(ηh,uh)
  
  m = DorflerMarking(0.8)
  I = Adaptivity.mark(m,η)
  
  method = Adaptivity.NVBRefinement(model)
  fmodel = Adaptivity.get_model(refine(method,model;cells_to_refine=I))

  error = sum(l2_norm(uh - u_exact,dΩ))
  return fmodel, uh, η, I, error
end

function test_amr(nsteps,order)
  model = LShapedModel(10)

  ϵ = 1e-2
  r(x) = ((x[1]-0.5)^2 + (x[2]-0.5)^2)^(1/2)
  u_exact(x) = 1.0 / (ϵ + r(x))

  vtk = false
  last_error = Inf
  for i in 1:nsteps
    fmodel, uh, η, I, error = amr_step(model,u_exact;order)
    if vtk
      is_refined = map(i -> ifelse(i ∈ I, 1, -1), 1:num_cells(model))
      Ω = Triangulation(model)
      writevtk(
        Ω,"tmp/model_$(i-1)",append=false,
        cellfields = [
          "uh" => uh,
          "η" => CellField(η,Ω),
          "is_refined" => CellField(is_refined,Ω),
          "u_exact" => CellField(u_exact,Ω),
        ],
      )
    end

    println("Error: $error, Error η: $(sum(η))")
    @test (i < 3) || (error < last_error)
    last_error = error
    model = fmodel
  end
end

############################################################################################

@testset "Dorfler marking" test_dorfler_marking()
@testset "AMR - Poisson" test_amr(20,2)

end # module