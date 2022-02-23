module Issue743

using Gridap
using LinearAlgebra
using Gridap.FESpaces
using Test

model = CartesianDiscreteModel((0,2π,0,2π),(10,10);isperiodic=(true,true))
order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe;constraint=:zeromean)
U = TrialFESpace(V)
degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
f(x) = cos(x[1])+sin(x[2])
a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
b(v) = ∫( v*f )*dΩ
op = AffineFEOperator(a,b,U,V)
uh = solve(op)
A = get_matrix(op)
@test cond(Matrix(A)) < 1.e5
@test abs(sum(∫(uh)dΩ)) < 1.e-10

#writevtk(Ω,"results",nsubcells=10,cellfields=["uh"=>uh])

end # module
