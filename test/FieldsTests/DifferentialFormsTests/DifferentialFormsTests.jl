module DifferentialFormsTests
# Lazy Field-level differential forms: DifferentialForm and its wrappers
# ExteriorDerivativeForm, CodifferentialForm, KoszulForm, PullbackForm,
# pushforward, and hodge_star_form. Single-point and multi-point evaluation
# through return_cache/evaluate!.

using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: DifferentialForm
using Gridap.Fields: ExteriorDerivativeForm
using Gridap.Fields: CodifferentialForm
using Gridap.Fields: KoszulForm
using Gridap.Fields: PullbackForm
using Gridap.Fields: hodge_star_form
using Test

fun(x) = 3.0*x[1]
f = GenericField(fun) # 3x
g = f*f               # 9x^2

D = 2
K = 1

ω1 = DifferentialForm{K,D}((f,g)) # 3x dx¹ + 9x² dx²
ω2 = DifferentialForm{K,D}((g,f)) # 9x² dx¹ + 3x dx²

# ── Single-point ──────────────────────────────────────────────────────────────

x = Point(1.0,1.0)

cache = return_cache(ω1,x)
ω1x = evaluate!(cache,ω1,x)
@test ω1x.data == (3.0,9.0)

cache = return_cache(ω2,x)
ω2x = evaluate!(cache,ω2,x)
@test ω2x.data == (9.0,3.0)

k = ω1 ∧ ω2
kx = evaluate(k,x)
@test kx == ω1x ∧ ω2x
@test kx.data == (-72.0,)

k = ω1 ∧ ω1
kx = evaluate(k,x)
@test kx == ω1x ∧ ω1x
@test kx.data == (0.0,)

# ── Multi-point: DifferentialForm ─────────────────────────────────────────────

xs = [Point(1.0,1.0), Point(2.0,1.0), Point(1.0,2.0)]

cache_v = return_cache(ω1, xs)
vals    = evaluate!(cache_v, ω1, xs)
@test vals isa AbstractVector
@test length(vals) == 3
@test vals[1].data == (3.0, 9.0)
@test vals[2].data == (6.0, 36.0)  # At (2,1): f = 3*2 = 6, g = f² = 36
@test vals[3].data == (3.0, 9.0)   # At (1,2): f = 3*1 = 3, g = 9

# cache reuse with different length
xs2   = [Point(0.5,0.5)]
vals2 = evaluate!(cache_v, ω1, xs2)
@test length(vals2) == 1
@test vals2[1].data == (1.5, 2.25)  # f=1.5, g=2.25

# ── Multi-point: ExteriorDerivativeForm ───────────────────────────────────────
# ω1 = 3x[1] dx¹ + 9x[1]² dx²   (components depend only on x[1])
# dω1 = (∂(9x[1]²)/∂x[1] - ∂(3x[1])/∂x[2]) dx¹∧dx² = (18x[1] - 0) dx¹∧dx²

dω1 = ExteriorDerivativeForm(ω1)

cache_d = return_cache(dω1, xs)
dvals   = evaluate!(cache_d, dω1, xs)
@test dvals isa AbstractVector
@test length(dvals) == 3
@test dvals[1].data[1] ≈ 18.0   atol=1e-12   # x[1]=1
@test dvals[2].data[1] ≈ 36.0   atol=1e-12   # x[1]=2
@test dvals[3].data[1] ≈ 18.0   atol=1e-12   # x[1]=1

# ── Multi-point: KoszulForm ───────────────────────────────────────────────────
# κ_x(ω) = ι_x ω.  At x=(1,1), ω1=(3,9): κ_x(ω) = 3*1 + 9*1 = 12

kf = KoszulForm(ω1)

cache_k = return_cache(kf, xs)
kvals   = evaluate!(cache_k, kf, xs)
@test kvals isa AbstractVector
@test length(kvals) == 3
# At (1,1): ι_{(1,1)}(3 dx¹ + 9 dx²) = 3*1 + 9*1 = 12
@test kvals[1].data[1] ≈ 12.0   atol=1e-12
# At (2,1): ω = (6,36), ι_{(2,1)} = 6*2 + 36*1 = 48
@test kvals[2].data[1] ≈ 48.0   atol=1e-12

# ── Single-point KoszulForm and koszul convenience constructor ────────────────

f1_fld = GenericField(x -> x[1])   # coefficient x¹
f2_fld = GenericField(x -> x[2])   # coefficient x²
ω_id   = DifferentialForm{1,2}((f1_fld, f2_fld))   # x¹ dx¹ + x² dx²
κω_fld = KoszulForm(ω_id)
@test typeof(κω_fld) <: KoszulForm{1,2}

x0 = Point(2.0, 3.0)
c = return_cache(κω_fld, x0)
κval = evaluate!(c, κω_fld, x0)
@test typeof(κval) <: DifferentialFormValue{0,2}
# κ(ω)(x₀) = x¹·x¹ + x²·x² = 4 + 9 = 13
@test κval.data[1] ≈ 13.0   atol=1e-12

κω2_fld = koszul(ω_id)
@test typeof(κω2_fld) <: KoszulForm{1,2}
c2 = return_cache(κω2_fld, x0)
@test evaluate!(c2, κω2_fld, x0).data[1] ≈ 13.0   atol=1e-12

# κ of a 2-form: κ(dx¹∧dx²)(x₀) = -x₀[2] dx¹ + x₀[1] dx²
f_vol = GenericField(_ -> 1.0)
ω_vol = DifferentialForm{2,2}((f_vol,))
κω_vol = KoszulForm(ω_vol)
c3 = return_cache(κω_vol, x0)
κval2 = evaluate!(c3, κω_vol, x0)
@test typeof(κval2) <: DifferentialFormValue{1,2}
# ι_{(2,3)}(dx¹∧dx²): d[1] = -x₀[2] = -3, d[2] = x₀[1] = 2
@test collect(κval2.data) ≈ [-3.0, 2.0]   atol=1e-12

# ── Field-level hodge_star_form ───────────────────────────────────────────────

fun1(x) = 3.0*x[1]*x[2]
fun2(x) = 5.0*x[2]^2
ω_fld  = DifferentialForm{1,2}((GenericField(fun1), GenericField(fun2)))  # 3xy dx¹ + 5y² dx²
star_ω = hodge_star_form(ω_fld)
@test typeof(star_ω) <: DifferentialForm{1,2}

x1 = Point(1.0, 2.0)
star_val     = evaluate(star_ω, x1)
expected_val = hodge_star(evaluate(ω_fld, x1))   # value-level star
@test collect(star_val.data) ≈ collect(expected_val.data)   atol=1e-12

# ── Field-level codifferential ────────────────────────────────────────────────

δω_fld = codifferential(ω_fld)
@test typeof(δω_fld) <: CodifferentialForm{1,2}

cδ   = return_cache(δω_fld, x1)
δval = evaluate!(cδ, δω_fld, x1)
@test typeof(δval) <: DifferentialFormValue{0,2}

# ω₁=3xy, ω₂=5y²; ⋆ω = (-ω₂, ω₁); d⋆ω = (∂ω₁/∂x + ∂ω₂/∂y) vol = (3y+10y) vol
# ⋆d⋆ω = 13y; δω = sign·⋆d⋆ω = -1·(13y).  At (1,2): -26
@test δval.data[1] ≈ -26.0   atol=1e-8

# ── Field-level pullback (PullbackForm) ───────────────────────────────────────

# Map φ: ℝ² → ℝ²  by φ(ξ) = (aξ¹, bξ²)
a, b = 2.0, 3.0
φ_fun(ξ) = VectorValue(a*ξ[1], b*ξ[2])
φ_fld    = GenericField(φ_fun)

# Form on the target: ω = x¹dx¹ + x²dx²
g1(x) = x[1]
g2(x) = x[2]
ω_target = DifferentialForm{1,2}((GenericField(g1), GenericField(g2)))

pb_fld = pullback(φ_fld, ω_target, Val(2))
@test typeof(pb_fld) <: PullbackForm{1,2,2}

ξ0  = Point(1.0, 1.0)
cpb = return_cache(pb_fld, ξ0)
pb_val = evaluate!(cpb, pb_fld, ξ0)
@test typeof(pb_val) <: DifferentialFormValue{1,2}
# φ*(x¹dx¹+x²dx²) = aξ¹·(adξ¹) + bξ²·(bdξ²) = a²ξ¹dξ¹ + b²ξ²dξ²
# At ξ=(1,1): (a², b²) = (4, 9)
@test collect(pb_val.data) ≈ [a^2, b^2]   atol=1e-8

# ── Field-level pushforward ───────────────────────────────────────────────────


function _pushforward(φ::Field, v::Field)
  Operation((J, u) -> J ⋅ u)(∇(φ), v)
end

v_fld = GenericField(ξ -> VectorValue(ξ[1], ξ[2]))   # v = ξ (identity vector field)
pf_fld = _pushforward(φ_fld, v_fld)
cpf    = return_cache(pf_fld, ξ0)
pf_val = evaluate!(cpf, pf_fld, ξ0)
@test typeof(pf_val) <: VectorValue{2}
# φ_*(ξ) = Jφ·ξ = diag(a,b)·(1,1) = (a,b)
@test collect(pf_val.data) ≈ [a, b]   atol=1e-8

end # module
