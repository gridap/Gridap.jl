using Gridap
using Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces, Gridap.Arrays
using Gridap.Fields

using BenchmarkTools

using Profile
using ProfileView

using FillArrays

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(40,40)))


order = 2
reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
reffe_RT = ReferenceFE(raviart_thomas,Float64,order-1)

V = FESpace(model,reffe)
V_RT = FESpace(model,reffe_RT)

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

a(u,v) = ∫(u⋅v)dΩ
b(u,v) = ∫((∇⋅u)*(∇⋅v))dΩ
c(u,v) = ∫(∇(u)⊙∇(v))dΩ

A1 = assemble_matrix(a,V,V)
A2 = assemble_matrix(a,V_RT,V_RT)

B1 = assemble_matrix(b,V,V)
B2 = assemble_matrix(b,V_RT,V_RT)

C1 = assemble_matrix(c,V,V)
C2 = assemble_matrix(c,V_RT,V_RT)

@benchmark assemble_matrix!(a,A1,V,V)
@benchmark assemble_matrix!(a,A2,V_RT,V_RT)
@benchmark assemble_matrix!(b,B1,V,V)
@benchmark assemble_matrix!(b,B2,V_RT,V_RT)
@benchmark assemble_matrix!(c,C1,V,V)
@benchmark assemble_matrix!(c,C2,V_RT,V_RT)

function p1()
  for i in 1:100
    assemble_matrix!(b,B2,V_RT,V_RT)
  end
end
ProfileView.@profview p1()


mon = V.fe_basis.cell_basis.values[1].fields
hdiv = V_RT.fe_basis.cell_basis.args[1].values[1].fields.qgrad



