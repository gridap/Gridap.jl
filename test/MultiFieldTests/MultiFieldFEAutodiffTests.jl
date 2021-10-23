module MultiFieldFEAutodiffTests

using Test

using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

Ω = Triangulation(model)
dΩ = Measure(Ω,2)

V1 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
V2 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
Y = MultiFieldFESpace([V1,V2])

r((uh,ph),(v,q)) = ∫( v*uh + uh*q + ph*q )dΩ
j(xh,dx,dy) = jacobian(xh->r(xh,dy),xh)

e((uh,ph)) = ∫( uh*uh + uh*ph + ph*ph )dΩ
g(xh,dy) = gradient(xh->e(xh),xh)
#h(xh,dx,dy) = hessian(xh->e(xh),xh)

dx = get_trial_fe_basis(Y)
dy = get_fe_basis(Y)
xh = FEFunction(Y,rand(num_free_dofs(Y)))

#display(r(xh,dy)[Ω][end])
#display(j(xh,dx,dy)[Ω][end])
#display(e(xh)[Ω][end])
#display(g(xh,dy)[Ω][end])
##display(h(xh,dx,dy)[Ω][end])

@test j(xh,dx,dy)[Ω][end][1,1] != nothing
@test j(xh,dx,dy)[Ω][end][2,1] != nothing
@test j(xh,dx,dy)[Ω][end][1,2] != nothing
@test j(xh,dx,dy)[Ω][end][2,2] != nothing
@test g(xh,dy)[Ω][end][1] != nothing
@test g(xh,dy)[Ω][end][2] != nothing

V1 = FESpace(model,ReferenceFE(lagrangian,Float64,2))
V2 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
Y = MultiFieldFESpace([V1,V2])

dx = get_trial_fe_basis(Y)
dy = get_fe_basis(Y)
xh = FEFunction(Y,rand(num_free_dofs(Y)))

@test j(xh,dx,dy)[Ω][end][1,1] != nothing
@test j(xh,dx,dy)[Ω][end][2,1] != nothing
@test j(xh,dx,dy)[Ω][end][1,2] != nothing
@test j(xh,dx,dy)[Ω][end][2,2] != nothing
@test g(xh,dy)[Ω][end][1] != nothing
@test g(xh,dy)[Ω][end][2] != nothing

eu(uh) = ∫( uh*uh )dΩ
ep(ph) = ∫( ph*ph )dΩ
eup((uh,ph)) = ∫( uh*uh + ph*ph )dΩ
geu(xh,dy) = gradient(xh->eu(xh),xh)
gep(xh,dy) = gradient(xh->ep(xh),xh)
geup(xh,dy) = gradient(xh->eup(xh),xh)

uh,ph=xh
geu_uh=geu(uh,dy)
gep_ph=gep(ph,dy)
geup_uh_ph=geup(xh,dy)

a=geu_uh[Ω]
b=gep_ph[Ω]
c=geup_uh_ph[Ω]

@test all(a .== map(x->x.array[1],c))
@test all(b .== map(x->x.array[2],c))


ru(uh,v) = ∫( v*uh*uh )dΩ
rp(ph,q) = ∫( q*ph*ph )dΩ
rup((uh,ph),(v,q)) = ∫( v*uh*uh + q*ph*ph )dΩ
ju(xh,dx,dy) = jacobian(xh->ru(xh,dy),xh)
jp(xh,dx,dy) = jacobian(xh->rp(xh,dy),xh)
jup(xh,dx,dy) = jacobian(xh->rup(xh,dy),xh)

uh=FEFunction(V1,get_free_dof_values(xh[1]))
ph=FEFunction(V2,get_free_dof_values(xh[2]))
du = get_trial_fe_basis(V1)
dv = get_fe_basis(V1)
dp = get_trial_fe_basis(V2)
dq = get_fe_basis(V2)


ju_uh=ju(uh,du,dv)
jp_ph=jp(ph,dp,dq)
jup_uh_ph=jup(xh,dx,dy)

a=ju_uh[Ω]
b=jp_ph[Ω]
c=jup_uh_ph[Ω]

@test all(a .== map(x->x.array[1,1],c))
@test all(b .== map(x->x.array[2,2],c))

end # module
