module DirichletFESpacesTests

using Test
using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.ReferenceFEs

n = 2
mesh = (n,n)
domain = (0,1,0,1)
order = 1
dirichlet_tags = [1,3,4,6,7]
model = CartesianDiscreteModel(domain, mesh)
Ω = Triangulation(model)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",dirichlet_tags)

reffe = ReferenceFE(lagrangian,Float64,order)
Vf = TestFESpace( model,reffe; conformity=:H1, labels=labels, dirichlet_tags="dirichlet")

Vd = DirichletFESpace(Vf)
cellmat = [rand(4,4) for cell in 1:num_cells(model)]
cellvec = [rand(4) for cell in 1:num_cells(model)]
cellmatvec = pair_arrays(cellmat,cellvec)
test_single_field_fe_space(Vd,cellmatvec,cellmat,cellvec,Ω)

Uf = TrialFESpace(Vf,1)

Ud = TrialFESpace(Vd)
test_single_field_fe_space(Ud,cellmatvec,cellmat,cellvec,Ω)

trian = Triangulation(model)
degree = 2*order
dΩ = Measure(trian,degree)

a(u,v) = ∫( u*v )*dΩ
l(v) = ∫( v*0 )*dΩ

op_ff = AffineFEOperator(a,l,Uf,Vf)
op_fd = AffineFEOperator(a,l,Ud,Vf)
op_df = AffineFEOperator(a,l,Uf,Vd)

xd = get_dirichlet_dof_values(Uf)

@test maximum(abs.(get_vector(op_ff) + get_matrix(op_fd)*xd)) < 1.0e-10

end # module
