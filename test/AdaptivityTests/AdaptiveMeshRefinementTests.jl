using Test
using Gridap, Gridap.Geometry, Gridap.Adaptivity
using DataStructures

# Marking tests

function test_marking()
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

test_marking()

# AMR tests

function LShapedModel(n)
  model = CartesianDiscreteModel((0,1,0,1),(n,n))
  cell_coords = map(mean,get_cell_coordinates(model))
  l_shape_filter(x) = (x[1] < 0.5) || (x[2] < 0.5)
  mask = map(l_shape_filter,cell_coords)
  return simplexify(DiscreteModelPortion(model,mask))
end

model = LShapedModel(10)

method = Adaptivity.NVBRefinement(model)
cells_to_refine = [collect(1:10)...,collect(20:30)...]
fmodel = refine(method,model;cells_to_refine)

writevtk(Triangulation(fmodel),"tmp/fmodel";append=false)

