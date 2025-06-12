module Issue614

using Gridap
using Gridap.CellData

domain = (0,2,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
const k = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},k)
reffe_p = ReferenceFE(lagrangian,Float64,k-1)
V₀ = TestFESpace(model, reffeᵤ, conformity=:H1)
U = TrialFESpace(V₀)
Q = TestFESpace(model, reffe_p, conformity=:C0)
P = TrialFESpace(Q)
Y = MultiFieldFESpace([V₀,Q])
X = MultiFieldFESpace([U,P])

Ω = Triangulation(model)

degree = 1
xh = interpolate_everywhere([VectorValue(0.0,1.0),x->x[1]],X)

Γc = BoundaryTriangulation(model)
dΓc = Measure(Γc, degree)
nΓc = get_normal_vector(Γc)

for i in 0:1000
    uh, ph = xh
    α = ph*nΓc
    CD, CL = ∑( ∫( α )dΓc )
end

end # module
