module FeecPoissonTests
# End-to-end validation of the differential-form bilinear forms for FEEC:
#  1. d_0form, d_1form at CellField level
#  2. Scalar Poisson:
#       a(u,v) = ∫(vol_coeff(d_0form(u) ∧ ⋆(d_0form(v)))) * dΩ
#     gives a machine-precision identical solution to ∫(∇u·∇v)*dΩ,
#     and the correct solution (L² error ≪ 1).

using Gridap
using Gridap.TensorValues
using Gridap.Geometry: CartesianDiscreteModel, Triangulation
using Gridap.CellData
using Gridap.Arrays: array_cache, getindex!
using Gridap.Fields
using Test

# ── 1. CellField level: d_0form, d_1form ─────────────────────────────────────

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0), (4,4))
Ω     = Triangulation(model)
pts   = Gridap.CellData.get_cell_points(Ω)

# d_0form on a scalar CellField: d(x¹+x²) should give (1,1) everywhere
u_scalar = CellField(x -> x[1] + x[2], Ω, PhysicalDomain())
du_cf = d_0form(u_scalar)
dux = du_cf(pts)
cdu = array_cache(dux)
duv = getindex!(cdu, dux, 1)
@test typeof(duv) <: AbstractArray{<:ExteriorFormValue{1,2}}
@test all(v -> collect(v.data) ≈ [1.0, 1.0], duv)

# d_1form on a VectorValue CellField: d(−x₂, x₁) = 2 dx¹∧dx² everywhere
u_vec = CellField(x -> VectorValue(-x[2], x[1]), Ω, PhysicalDomain())
ddu_cf = d_1form(u_vec)
ddux = ddu_cf(pts)
cddu = array_cache(ddux)
dduv = getindex!(cddu, ddux, 1)
@test typeof(dduv) <: AbstractArray{<:ExteriorFormValue{2,2}}
@test all(v -> vol_coeff(v) ≈ 2.0, dduv)

# ── 2. Scalar Poisson: differential-form weak form vs standard Gridap ────────
# −Δu = f on [0,1]²,  u = 0 on ∂Ω
# Exact solution:  u = sin(πx)sin(πy),  f = 2π²sin(πx)sin(πy)

model8 = CartesianDiscreteModel((0.0,1.0,0.0,1.0), (8,8))
Ω8     = Triangulation(model8)
dΩ8    = Measure(Ω8, 4)

V = FESpace(model8, ReferenceFE(QUAD, lagrangian, Float64, 2);
            conformity=:H1, dirichlet_tags="boundary")
U = TrialFESpace(V, 0.0)

f_cf = CellField(x -> 2π^2*sin(π*x[1])*sin(π*x[2]), Ω8, PhysicalDomain())
l(v) = ∫( f_cf * v ) * dΩ8

# Differential-form bilinear form: a(u,v) = ∫(d₀u ∧ ⋆d₀v)
a_form(u,v) = ∫( vol_coeff(d_0form(u) ∧ hodge_star(d_0form(v))) ) * dΩ8

# Standard Gridap bilinear form for comparison
a_std(u,v) = ∫( ∇(u) ⋅ ∇(v) ) * dΩ8

op_form = AffineFEOperator(a_form, l, U, V)
op_std  = AffineFEOperator(a_std,  l, U, V)

uh_form = solve(op_form)
uh_std  = solve(op_std)

# Both formulations must give machine-precision identical solutions
diff = uh_form - uh_std
@test sqrt(sum(∫(diff * diff) * dΩ8)) < 1e-12

# Solution must converge to the exact answer (order-2, h=1/8 → expect ≈ 4e-4)
u_ex = CellField(x -> sin(π*x[1])*sin(π*x[2]), Ω8, PhysicalDomain())
e    = uh_form - u_ex
@test sqrt(sum(∫(e * e) * dΩ8)) < 5e-4

end # module
