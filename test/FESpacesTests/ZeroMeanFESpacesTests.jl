module ZeroMeanFESpacesTests

using Test
using FillArrays
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.FESpaces
using Gridap.CellData
using Gridap.ReferenceFEs

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 2

trian = get_triangulation(model)
degree = order
dΩ = Measure(trian,degree)

_V = FESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:L2)

V = ZeroMeanFESpace(_V,dΩ)

cellmat = [rand(4,4) for cell in 1:num_cells(model)]
cellvec = [rand(4) for cell in 1:num_cells(model)]
cellmatvec = pair_arrays(cellmat,cellvec)
test_single_field_fe_space(V,cellmatvec,cellmat,cellvec,trian)

U = TrialFESpace(V)
test_single_field_fe_space(U,cellmatvec,cellmat,cellvec,trian)
@test isa(U,ZeroMeanFESpace)

# Interpolate a function with non-zero mean
f(x) = sin(4*pi*(x[1]+x[2]^2)) + 3
uh = interpolate(f, U)
@test abs(sum(∫(uh)*dΩ)) < 1.0e-10

# Interpolate a function with zero mean
ĝ(x) = x[1]^2 + x[2]
g_mean = sum(∫(ĝ)*dΩ)/sum(∫(1)*dΩ)
g(x) = ĝ(x) - g_mean
vh = interpolate(g, U)
eh = vh - g
@test abs(sum(∫(vh)*dΩ)) < 1.0e-10
@test sum(∫(eh*eh)*dΩ) < 1.0e-10

# Check Stokes/Navier-Stokes problem properties

u_ex(x) = VectorValue(x[2],-x[1])
p_ex(x) = x[1] + 2*x[2]
p_mean = sum(∫(p_ex)*dΩ)

order = 2
reffe_u = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
reffe_p = ReferenceFE(lagrangian,Float64,order-1;space=:P)

V  = FESpace(model,reffe_u;dirichlet_tags="boundary")
U  = TrialFESpace(V,u_ex)
Q  = FESpace(model,reffe_p;conformity=:L2)
Q0 = FESpace(model,reffe_p;conformity=:L2,constraint=:zeromean)
@test isa(Q0,ZeroMeanFESpace)

X = MultiFieldFESpace([U,Q0])
Y = MultiFieldFESpace([V,Q0])

ph_i = interpolate(p_ex, Q0)
@test abs(sum(∫(ph_i)*dΩ)) < 1.0e-10
@test abs(sum(∫(ph_i - p_ex + p_mean)*dΩ)) < 1.0e-10

a((u,p),(v,q)) = ∫(∇(u)⊙∇(v) - (∇⋅v)*p - q*(∇⋅u))dΩ
l((v,q)) = a((u_ex,p_ex),(v,q))

op = AffineFEOperator(a,l,X,Y)
uh, ph = solve(op)

l2_error(u,v) = sqrt(sum(∫((u-v)⋅(u-v))*dΩ))
@test l2_error(uh,u_ex) < 1.0e-10
@test abs(sum(∫(ph)*dΩ)) < 1.0e-10
@test l2_error(ph,ph_i) < 1.0e-10

b(u,q) = ∫(q*(∇⋅u))dΩ
B = assemble_vector(q -> b(uh,q),Q)
B0 = assemble_vector(q -> b(uh,q),Q0)
@test abs(sum(B)) < 1.0e-10
@test abs(sum(B0)) < 1.0e-10

uh_ex = interpolate(u_ex, U)
@test l2_error(uh,uh_ex) < 1.0e-10
conv(u,∇u) = (∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)
c(u,v) = ∫(v⊙(conv∘(u,∇(u))))dΩ
dc(u,du,dv) = ∫(dv⊙(dconv∘(du,∇(du),u,∇(u))))dΩ
jac((u,p),(du,dp),(dv,dq)) = a((du,dp),(dv,dq)) + dc(u,du,dv)
res((u,p),(dv,dq)) = a((u,p),(dv,dq)) + c(u,dv) - a((u_ex,p_ex),(dv,dq)) - c(uh_ex,dv)

op = FEOperator(res,jac,X,Y)
uh, ph = solve(op)

@test l2_error(uh,u_ex) < 1.0e-10
@test abs(sum(∫(ph)*dΩ)) < 1.0e-10
@test l2_error(ph,ph_i) < 1.0e-10

B = assemble_vector(q -> b(uh,q),Q)
B0 = assemble_vector(q -> b(uh,q),Q0)
@test abs(sum(B)) < 1.0e-10
@test abs(sum(B0)) < 1.0e-10

# Moment based bases
#
# `modalC0` is not a partition of unity from order 2 on, so the DOF values of 𝔥
# are not all ones: the volume is not the sum of the `∫φᵢ`, and the DOF the
# space fixes cannot be the last one, on which they vanish.

# A nodal basis keeps the single scalar shift
@test ZeroMeanFESpace(_V,dΩ).h_free isa Ones

# A space with Dirichlet DOFs is rejected
@test_throws ArgumentError ZeroMeanFESpace(
  FESpace(model,ReferenceFE(lagrangian,Float64,1);dirichlet_tags="boundary"),dΩ)

# A vector valued space is rejected
@test_throws "scalar valued spaces only" ZeroMeanFESpace(
  FESpace(model,ReferenceFE(lagrangian,VectorValue{2,Float64},1)),dΩ)

for p_order in (2,3)
  _Vm = FESpace(model,ReferenceFE(QUAD,modalC0,Float64,p_order);conformity=:L2)
  Vm  = ZeroMeanFESpace(_Vm,dΩ)
  @test !isa(Vm.h_free,Ones)
  @test Vm.vol ≈ sum(∫(1)*dΩ)
  # the DOF values of 𝔥 vanish on the last DOF, the one a partition of unity
  # fixes, so the space has to fix another one
  @test abs(get_free_dof_values(interpolate(1.0,_Vm))[end]) < 1.0e-12
  @test Vm.h_diri != 0

  vhm = interpolate(g,TrialFESpace(Vm))
  ehm = vhm - g
  @test abs(sum(∫(vhm)*dΩ)) < 1.0e-10
  @test sum(∫(ehm*ehm)*dΩ) < 1.0e-10
end

# The invariant the construction checks: the vector of ones does not represent
# the mean functional of a moment based basis, the interpolation of 1 does.
_Vm = FESpace(model,ReferenceFE(QUAD,modalC0,Float64,2);conformity=:L2)
int_φi_m = assemble_vector(v->∫(v)*dΩ,_Vm)
@test !FESpaces._represents_mean_functional(_Vm,dΩ,int_φi_m,
        FEFunction(_Vm,ones(num_free_dofs(_Vm))))
@test FESpaces._represents_mean_functional(_Vm,dΩ,int_φi_m,interpolate(1.0,_Vm))

# Pure Neumann Poisson, whose solution the zero mean value selects. `modalC0`
# of order 2 spans the same space as `lagrangian` of order 2, so both bases
# give the same solution.
u_N(x) = cos(pi*x[1])
f_N(x) = pi^2*cos(pi*x[1])
dΩ_N = Measure(trian,8)

function zero_mean_poisson(reffe)
  V0 = FESpace(model,reffe;conformity=:H1,constraint=:zeromean)
  op = AffineFEOperator((u,v) -> ∫( ∇(u)⋅∇(v) )*dΩ_N,
                        v -> ∫( f_N*v )*dΩ_N, TrialFESpace(V0), V0)
  solve(op)
end

uh_nodal  = zero_mean_poisson(ReferenceFE(lagrangian,Float64,2))
uh_moment = zero_mean_poisson(ReferenceFE(QUAD,modalC0,Float64,2))
@test abs(sum(∫(uh_nodal)*dΩ_N))  < 1.0e-10
@test abs(sum(∫(uh_moment)*dΩ_N)) < 1.0e-10
e_N = uh_nodal - uh_moment
@test sqrt(sum(∫( e_N*e_N )*dΩ_N)) < 1.0e-10
e_ex = uh_nodal - u_N
@test sqrt(sum(∫( e_ex*e_ex )*dΩ_N)) < 1.0e-2

end # module
