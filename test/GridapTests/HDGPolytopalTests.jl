module HDGPolytopalTests

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField, Gridap.Polynomials
using Gridap.CellData, Gridap.Fields, Gridap.ReferenceFEs, Gridap.Helpers
using Gridap.Arrays

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

#################
# HDG+ variant
#################
#u(x) = sin(2*π*x[1])*sin(2*π*x[2])
u(x) = x[1] + x[2]
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

D = 3
nc = Tuple(fill(4, D))
domain = Tuple(repeat([0, 1], D))
model = simplexify(CartesianDiscreteModel(domain,nc);positive=true)

if (D == 2)
  vmodel = Gridap.Geometry.voronoi(model)
else
  vmodel = Gridap.Geometry.PolytopalDiscreteModel(model)
end

# writevtk(vmodel,"tmp/polygonal_model")

Ω = Triangulation(ReferenceFE{D}, vmodel)
Γ = Triangulation(ReferenceFE{D-1}, vmodel)

ptopo = Geometry.PatchTopology(vmodel)
Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

# Reference FEs
order = 1
V = FESpaces.PolytopalFESpace(Ω, VectorValue{D, Float64}, order; space=:P)
Q = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P)

if D == 2
  reffeM = ReferenceFE(lagrangian, Float64, order; space=:P) 
  M = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")
else
  M = FESpaces.PolytopalFESpace(Γ, Float64, order; space=:P, dirichlet_tags="boundary")
end
N = TrialFESpace(M, u)

mfs = MultiField.BlockMultiFieldStyle(2,(2,1))
X = MultiFieldFESpace([V, Q, N];style=mfs)

degree = 2*(order+1)
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)

τT = CellField(1 ./ get_array(∫(1)dΓp) .^(D-1), Γp) # HDG stab parameter
n = get_normal_vector(Γp)
Πn(u) = u⋅n
Π(u) = change_domain(u,Γp,DomainStyle(u))

PΓ = projection_operator(M, Γp, dΓp)

Pqhn(qh,uh,sh) = Πn(qh) + τT * (PΓ(uh)-sh) # double vector valued function on mesh / Transmission condition
a((qh,uh,sh),(vh,wh,lh)) = ∫( qh⋅vh - uh*(∇⋅vh) - qh⋅∇(wh) )dΩp + ∫(sh*Πn(vh))dΓp +
                           ∫(Pqhn(qh,uh,sh)*(Π(wh) + lh))dΓp 
l((vh,wh,lh)) = ∫( f*wh )*dΩp

op = MultiField.StaticCondensationOperator(ptopo,X,a,l)
qh, uh, sh = solve(op)

dΩ = Measure(Ω,degree)
eh = uh - u
l2_uh = sqrt(sum(∫(eh⋅eh)*dΩ))
@test l2_uh < 1e-10

end # module