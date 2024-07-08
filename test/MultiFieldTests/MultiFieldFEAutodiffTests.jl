module MultiFieldFEAutodiffTests

using Test

using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

using ForwardDiff

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
h(xh,dx,dy) = hessian(xh->e(xh),xh)

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
@test h(xh,dx,dy)[Ω][end][1,1] != nothing
@test h(xh,dx,dy)[Ω][end][2,1] != nothing
@test h(xh,dx,dy)[Ω][end][1,2] != nothing
@test h(xh,dx,dy)[Ω][end][2,2] != nothing

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
@test h(xh,dx,dy)[Ω][end][1,1] != nothing
@test h(xh,dx,dy)[Ω][end][2,1] != nothing
@test h(xh,dx,dy)[Ω][end][1,2] != nothing
@test h(xh,dx,dy)[Ω][end][2,2] != nothing

eu(uh) = ∫( uh*uh )dΩ
ep(ph) = ∫( ph*ph )dΩ
ez((uh,ph)) = ∫( 0*uh*ph )dΩ
eup((uh,ph)) = ∫( uh*uh + ph*ph )dΩ
geu(xh,dy) = gradient(xh->eu(xh),xh)
gep(xh,dy) = gradient(xh->ep(xh),xh)
geup(xh,dy) = gradient(xh->eup(xh),xh)
heu(xh,dx,dy) = hessian(xh->eu(xh),xh)
hep(xh,dx,dy) = hessian(xh->ep(xh),xh)
hez(xh,dx,dy) = hessian(xh->ez(xh),xh)
heup(xh,dx,dy) = hessian(xh->eup(xh),xh)

uh,ph=xh
geu_uh=geu(uh,dy)
gep_ph=gep(ph,dy)
geup_uh_ph=geup(xh,dy)

a=geu_uh[Ω]
b=gep_ph[Ω]
c=geup_uh_ph[Ω]

@test all(a .== map(x->x.array[1],c))
@test all(b .== map(x->x.array[2],c))

heu_uh=heu(uh,dy,dy)
hep_ph=hep(ph,dy,dy)
hez_uh=hez(xh,dy,dy)
heup_uh_ph=heup(xh,dy,dy)

a=heu_uh[Ω]
b=hep_ph[Ω]
c=hez_uh[Ω]
d=heup_uh_ph[Ω]

@test all(a .== map(x->x.array[1,1],d))
@test all(b .== map(x->x.array[2,2],d))
@test all(map(x->x.array[2,1],c) .== map(x->x.array[2,1],d))
@test all(map(x->x.array[1,2],c) .== map(x->x.array[1,2],d))

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

# AD tests for Multifiled functionals involving SkeletonTriangulation
# compared with direct ForwardDiff results

reffeV = ReferenceFE(lagrangian,Float64,2)
reffeQ = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(model,reffeV,conformity=:L2)
Q = TestFESpace(model,reffeQ,conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

xh = FEFunction(X,rand(num_free_dofs(X)))

Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ,2)
n_Λ = get_normal_vector(Λ)

g_Λ((uh,ph)) = ∫( mean(uh) + mean(ph) + mean(uh)*mean(ph) )dΛ
f_Λ((uh,ph)) = ∫( mean(uh*uh) + mean(uh*ph) + mean(ph*ph) )dΛ
a_Λ((uh,ph)) = ∫( - jump(uh*n_Λ)⊙mean(∇(ph))
                  - mean(∇(uh))⊙jump(ph*n_Λ)
                  + jump(uh*n_Λ)⊙jump(ph*n_Λ) )dΛ

function f_uh_free_dofs(f,xh,θ)
  dir_u = similar(xh[1].dirichlet_values,eltype(θ))
  dir_p = similar(xh[2].dirichlet_values,eltype(θ))
  X = xh.fe_space
  θ_uh = restrict_to_field(X,θ,1)
  θ_ph = restrict_to_field(X,θ,2)
  uh = FEFunction(X[1],θ_uh,dir_u)
  ph = FEFunction(X[2],θ_ph,dir_p)
  xh = MultiFieldFEFunction(θ,xh.fe_space,[uh,ph])
  sum(f(xh))
end
# can also do the above by constructing a MultiFieldCellField

g_Λ_(θ) = f_uh_free_dofs(g_Λ,xh,θ)
f_Λ_(θ) = f_uh_free_dofs(f_Λ,xh,θ)
a_Λ_(θ) = f_uh_free_dofs(a_Λ,xh,θ)

θ = get_free_dof_values(xh)

# check if the evaluations are working
@test sum(f_Λ(xh)) == f_Λ_(θ)
@test sum(a_Λ(xh)) == a_Λ_(θ)
@test sum(g_Λ(xh)) == g_Λ_(θ)

f_Ω((uh,ph)) = ∫(uh*uh + ph*uh + uh*ph)*dΩ
f_Ω_(θ) = f_uh_free_dofs(f_Ω,xh,θ)
gridapgradf_Ω = assemble_vector(gradient(f_Ω,xh),X)
fdgradf_Ω = ForwardDiff.gradient(f_Ω_,θ)
test_array(gridapgradf_Ω,fdgradf_Ω,≈)

gridapgradf = assemble_vector(gradient(f_Λ,xh),X)
fdgradf = ForwardDiff.gradient(f_Λ_,θ)
test_array(gridapgradf,fdgradf,≈)

gridapgradg = assemble_vector(gradient(g_Λ,xh),X)
fdgradg = ForwardDiff.gradient(g_Λ_,θ)
test_array(gridapgradg,fdgradg,≈)

gridapgrada = assemble_vector(gradient(a_Λ,xh),X)
fdgrada = ForwardDiff.gradient(a_Λ_,θ)
test_array(gridapgrada,fdgrada,≈)

end # module
