module Issue770

using Gridap
using Test

domain = (0,1,0,1)
cells = (2,2)
model = CartesianDiscreteModel(domain,cells) |> simplexify
Ω = Interior(model)
dΩ = Measure(Ω,4)
m = 2
u(x) = VectorValue(sum(x)^m, sum(x)^m)
reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},m)
V = FESpace(Ω,reffe)
uh = interpolate(u,V)
t = Δ(u) - Δ(uh)
@test sum(∫(t⋅t)dΩ) < 1.0e-12

s(x) = sum(x)^m
reffe = ReferenceFE(lagrangian,Float64,m)
V = FESpace(Ω,reffe)
sh = interpolate(s,V)
t = Δ(s) - Δ(sh)
@test sum(∫(t⋅t)dΩ) < 1.0e-12

end # module
