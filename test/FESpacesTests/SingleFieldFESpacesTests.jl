module SingleFieldFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
reffe = ReferenceFE(:Lagrangian,Float64,order)
V0 = FESpace(model,reffe,dirichlet_tags=["tag_24","tag_25"])

matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(V0,matvecdata,matdata,vecdata)

f(x) = sin(4*pi*(x[1]-x[2]^2))+1

fh = interpolate_everywhere(f, V0)

fh = interpolate_dirichlet(f, V0)

dirichlet_values = compute_dirichlet_values_for_tags(V0,[1,2])

r=[2,2,1,1,2,2,1,1,2,2,2,2,1,2,1,1,1,1,2,2,2,2,1,2,1,1,1,1,2,2,1,1,2,2,1,2,1,1,2,2,1,2,1,1,2,2,1,2,1,1]

@test dirichlet_values == r

dirichlet_values = compute_dirichlet_values_for_tags(V0,2)
@test dirichlet_values == fill(2,num_dirichlet_dofs(V0))

free_values = zero_free_values(V0)
uh = FEFunction(V0,free_values,dirichlet_values)

# With complex dofs values
V0 = FESpace(model,reffe,dirichlet_tags=["tag_24","tag_25"],vector_type=Vector{ComplexF64})
@test get_vector_type(V0) == Vector{ComplexF64}
f(x) = x[1] + x[2]*im
fh = interpolate(f, V0)
fh = interpolate_dirichlet(f, V0)
fh = interpolate_everywhere(f, V0)

@test isa(real(fh),FEFunction)
@test isa(imag(fh),FEFunction)

f_real(x) = real(f(x))
f_imag(x) = imag(f(x))

Ω = Triangulation(model)
dΩ = LebesgueMeasure(Ω,2)

tol = 1e-9
@test sqrt(sum(∫( abs2(f - fh) )*dΩ)) < tol
@test sqrt(sum(∫( abs2(f_imag - imag(fh)) )*dΩ)) < tol
@test sqrt(sum(∫( abs2(f_real - real(fh)) )*dΩ)) < tol
@test sqrt(sum(∫( abs2(real(f - fh)) )*dΩ)) < tol
@test sqrt(sum(∫( abs2(imag(f - fh)) )*dΩ)) < tol

zh = zero(V0)
@test isa(get_free_values(zh),Vector{ComplexF64})


#using Gridap.Visualization
#trian = get_triangulation(model)
#writevtk(trian,"trian",cellfields=["fh_real"=>real(fh),"fh_imag"=>imag(fh)])



#using Gridap.Visualization
#
#trian = get_triangulation(model)
#
#writevtk(model,"model")
#
#writevtk(trian,"trian",cellfields=["fh" => fh, "uh" => uh])

#
#
#strian = SkeletonTriangulation(model)
#
#fh_Γ = restrict(get_array(fh),strian)
#
#writevtk(strian,"strian",nsubcells=10,cellfields=["jump_fh" => jump(fh_Γ)])

end # module
