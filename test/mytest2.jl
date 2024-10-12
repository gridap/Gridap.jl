
using SparseArrays, LinearAlgebra
using Gridap, Gridap.FESpaces, Gridap.Algebra, Gridap.Adaptivity

using GridapDistributed, PartitionedArrays

sol(x) = x

np = (2,2)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

cmodel = CartesianDiscreteModel(ranks,np,(0,1,0,1),(4,4))
fmodel = refine(cmodel)

reffe_u = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
reffe_j = ReferenceFE(raviart_thomas, Float64, 0)

Vh = TestFESpace(fmodel, reffe_u, dirichlet_tags="boundary")
Dh = TestFESpace(fmodel, reffe_j, dirichlet_tags="boundary")
Uh = TrialFESpace(Vh, sol)
Jh = TrialFESpace(Dh, sol)

VH = TestFESpace(cmodel, reffe_u, dirichlet_tags="boundary")
DH = TestFESpace(cmodel, reffe_j, dirichlet_tags="boundary")
UH = TrialFESpace(VH, sol)
JH = TrialFESpace(DH, sol)

f(x) = VectorValue(x[1]-x[2], x[1]+x[2])

Xh,XH,Yh,YH = Uh,UH,Vh,VH
#Xh,XH,Yh,YH = Jh,JH,Dh,DH

fv_h = zero_free_values(Xh)
fv_H = zero_free_values(XH)

dv_h = zero_dirichlet_values(Xh)
dv_H = zero_dirichlet_values(XH)

fill!(fv_H,1.0)
uH = FEFunction(XH,fv_H,dv_H)
uh = interpolate!(uH,fv_h,Xh)

fv_h_bis = zero_free_values(Yh)
fv_H_bis = zero_free_values(YH)

fill!(fv_H_bis,1.0)
vH = FEFunction(XH,fv_H_bis,dv_H)
vh = interpolate!(uH,fv_h_bis,Xh)

fv_H_bis == fv_H
