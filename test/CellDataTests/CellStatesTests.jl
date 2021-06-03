module CellStatesTests

using Test
using FillArrays
using Gridap.Helpers
using Gridap.Fields
using Gridap.CellData
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry

domain = (0,1,0,1)
cells = (3,3)
model = CartesianDiscreteModel(domain,cells)

Ω = Triangulation(model)
dΩ = Measure(Ω,2)

f = CellState{Float64}(undef,dΩ)
@test f.values[1] !== f.values[2]
updater(f) = true, 1.0
update_state!(updater,f)
@test sum(∫(f)*dΩ) ≈ 1

f = CellState(1,dΩ)
@test f.values[1] !== f.values[2]
@test sum(∫(f)*dΩ) ≈ 1

f = CellState(2.0,dΩ)
g = CellState(VectorValue(4.0,4.0),dΩ)
h = CellState(2.0,dΩ)
updater(f,g,h) = true, 0.5*f, 0.25*g
update_state!(updater,f,g,h)
@test sum(∫(f)*dΩ) ≈ 1
@test sum(∫(g)*dΩ) ≈ VectorValue(1,1)
@test sum(∫(h)*dΩ) ≈ 2

f = CellState(2.0,dΩ)
g = CellState(VectorValue(4.0,4.0),dΩ)
h = CellField(2.0,Ω)
updater(f,g,h) = true, 0.5*f, 0.25*g
update_state!(updater,f,g,h)
@test sum(∫(f)*dΩ) ≈ 1
@test sum(∫(g)*dΩ) ≈ VectorValue(1,1)
@test sum(∫(h)*dΩ) ≈ 2

r1 = f*h
@test sum(∫(r1)*dΩ) ≈ 2

σ(f,g,h) = f*g*h
r2 = σ∘(f,g,h)
@test sum(∫(r2)*dΩ) ≈ VectorValue(2,2)

r = r1*r2*h
@test sum(∫(r)*dΩ) ≈ VectorValue(8,8)

#x = get_cell_points(dΩ)
#using Gridap.Visualization
#writevtk(x,"x",cellfields=["f"=>f,"g"=>g])

end # module
