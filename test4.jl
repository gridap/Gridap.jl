using Gridap
using FillArrays
using Test 

f(x)=x[1]^2+x[2]

modelH=CartesianDiscreteModel((0,1,0,1),(1,1))
modelh=Gridap.Adaptivity.refine(modelH,2)

order  = 2
reffeH = ReferenceFE(lagrangian,Float64,order)
VH     = TestFESpace(modelH,reffeH; conformity=:H1)#,dirichlet_tags="boundary")
UH     = TrialFESpace(VH)

duH=get_trial_fe_basis(UH)

Ωh   = Triangulation(modelh)
dΩh  = Measure(Ωh,2*order)
ΩH   = Triangulation(modelH)
dΩH  = Measure(ΩH,2*order)

Vh=Gridap.LinearizedFESpace(modelH,reffeH;
                            conformity=:H1)#,
                            #dirichlet_tags="boundary")
                                    
Uh=TrialFESpace(Vh,0.0)

Ωh   = Triangulation(modelh)
dΩh  = Measure(Ωh,2*order)
ΩH   = Triangulation(modelH)
dΩH  = Measure(ΩH,2*order)

dΩHh = Measure(ΩH,Ωh,2*order)

adΩHh(UH,vh)=∫(UH*vh)dΩHh
AdΩHh=assemble_matrix(adΩHh,UH,Vh)

adΩH(UH,vh)=∫(UH*vh)dΩH
AdΩH=assemble_matrix(adΩH,UH,Vh)

norm(AdΩH-AdΩHh) # High error with the "wrong" measure

l(vh)=∫(1.0*vh)dΩHh
b1ΩHh=assemble_vector(l,Vh)

l(vh)=∫(1.0*vh)dΩH
b1ΩH=assemble_vector(l,Vh)

norm(b1ΩHh-b1ΩH) # High error with the "wrong" measure

reffe = ReferenceFE(lagrangian,Float64,1)
Vhl   = TestFESpace(modelh,reffe; conformity=:H1)#,dirichlet_tags="boundary")
l(vh)=∫(1.0*vh)dΩh
b2=assemble_vector(l,Vhl)

norm(b1ΩHh-b2) # Eps error in the "right" measure
norm(b1ΩH-b2)  # High error with the "wrong" measure

# Test L2-projection
adΩHh(UH,vh)=∫(UH*vh)dΩHh
lΩHh(vh)    =∫(f*vh)dΩHh
op=AffineFEOperator(adΩHh,lΩHh,UH,Vh)
fh=solve(op)
f(x)=x[1]^2+x[2]
eh=sum(∫((f-fh)*(f-fh))dΩHh)
@test eh < 1.0e-12

uH=get_trial_fe_basis(UH)
∇(uH)

vH=get_fe_basis(VH)
∇(vH)


vh=get_fe_basis(Vh)
grad_vh=∇(vh)
grad_vh_data=Gridap.CellData.get_data(grad_vh)
lazy_map(evaluate,grad_vh_data,Fill([Point(0.0,0.0)],4))

vhl=get_fe_basis(Vhl)
grad_vhl=∇(vhl)
grad_vhl_data=Gridap.CellData.get_data(grad_vhl)
lazy_map(evaluate,grad_vhl_data,Fill([Point(0.0,0.0)],4))


adΩHh(UH,vh)=∫(∇(UH)⋅∇(vh))dΩHh
lΩHh(vh)    =∫(f*vh)dΩHh
op=AffineFEOperator(adΩHh,lΩHh,UH,Vh)

