module DirichletFESpacesTests

using Test
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces

n = 2
mesh = (n,n)
domain = (0,1,0,1)
order = 1
dirichlet_tags = [1,3,4,6,7]
model = CartesianDiscreteModel(domain, mesh)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",dirichlet_tags)

Vf = TestFESpace(
  reffe=:Lagrangian, order=order, valuetype=Float64,
  conformity=:H1, model=model, labels=labels, dirichlet_tags="dirichlet")

Vd = DirichletFESpace(Vf)
matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(Vd,matvecdata,matdata,vecdata)

Uf = TrialFESpace(Vf,1)

Ud = TrialFESpace(Vd)
test_single_field_fe_space(Ud,matvecdata,matdata,vecdata)

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

a(v,u) = v*u
l(v) = v*0
t_Ω = AffineFETerm(a,l,trian,quad)

op_ff = AffineFEOperator(Uf,Vf,t_Ω)
op_fd = AffineFEOperator(Ud,Vf,t_Ω)
op_df = AffineFEOperator(Uf,Vd,t_Ω)

xd = get_dirichlet_values(Uf)

@test maximum(abs.(get_vector(op_ff) + get_matrix(op_fd)*xd)) < 1.0e-10

end # module
