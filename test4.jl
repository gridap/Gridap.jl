using Gridap
using FillArrays

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

duH=get_trial_fe_basis(UH)
dvh=get_fe_basis(Vh)

dvh_data=Gridap.CellData.get_data(dvh)
lazy_map(evaluate,dvh_data,Fill([Point(0.0,0.75)],1))[1]

Ωh   = Triangulation(modelh)
dΩh  = Measure(Ωh,2*order)
ΩH   = Triangulation(modelH)
dΩH  = Measure(ΩH,2*order)

xH=Gridap.CellData.get_cell_points(dΩH.quad)
dΩHh=Measure(ΩH,Ωh,2*order)

adΩHh(UH,vh)=∫(UH*vh)dΩHh
AdΩHh=assemble_matrix(a,UH,Vh)

adΩH(UH,vh)=∫(UH*vh)dΩH
AdΩH=assemble_matrix(a,UH,Vh)

AdΩH-AdΩHh

l(vh)=∫(1.0*vh)dΩHh
b1ΩHh=assemble_vector(l,Vh)

l(vh)=∫(1.0*vh)dΩH
b1ΩH=assemble_vector(l,Vh)

norm(b1ΩHh-b1ΩH) # High error with the "wrong" measure

reffe = ReferenceFE(lagrangian,Float64,1)
Vhl   = TestFESpace(modelh,reffe; conformity=:H1)#,dirichlet_tags="boundary")
l(vh)=∫(1.0*vh)dΩh
b2=assemble_vector(l,Vhl)

norm(b1ΩHh-b2)
norm(b1ΩH-b2) # High error with the "wrong" measure

using LinearAlgebra
eigvals(Array(A))


