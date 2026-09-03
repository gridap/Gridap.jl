module GridapSymbolicsTests
# Tests for the GridapSymbolicsExt extension: symbolic coordinates, symbolic
# exterior derivative/codifferential/Lie derivative on Symbolics.Num-valued
# DifferentialFormValues, the Koszul homotopy identity, and the print_forms
# symbolic display of the rotating/trimmed PΛ bases.

using Symbolics
using Gridap
using Gridap.TensorValues
using Gridap.Polynomials
using Test

# ── Symbolic coordinates ─────────────────────────────────────────────────────

x¹, x² = symbolic_coordinates(2)
@test x¹ isa Symbolics.Num
@test length(symbolic_coordinates(3)) == 3

v2 = VectorValue(x¹, x²)

# ── Symbolic scalar multiplication ───────────────────────────────────────────

ωs = DifferentialFormValue{1,2}((1.0, 2.0))
sω = x¹ * ωs
@test typeof(sω) <: DifferentialFormValue{1,2,Symbolics.Num}
@test isequal(Symbolics.simplify(sω.data[2] - 2x¹), Symbolics.simplify(Symbolics.Num(0)))
@test isequal((ωs * x¹).data[1], sω.data[1])

# ── Symbolic exterior derivative (Cartesian) ─────────────────────────────────

ω0 = DifferentialFormValue{0,2}((x¹*x²,))
dω0 = exterior_derivative(ω0)
@test typeof(dω0) <: DifferentialFormValue{1,2}
@test isequal(Symbolics.simplify(dω0.data[1] - x²), Symbolics.simplify(Symbolics.Num(0)))
@test isequal(Symbolics.simplify(dω0.data[2] - x¹), Symbolics.simplify(Symbolics.Num(0)))

# ── Exterior derivative with explicit vars (barycentric λ) and d² = 0 ────────

λ₁, λ₂, λ₃ = Symbolics.variables(:λ, 1:3)
λ = (λ₁, λ₂, λ₃)

ωb = DifferentialFormValue{0,3}((λ₁*λ₂,))
dωb = exterior_derivative(ωb, λ)
@test typeof(dωb) <: DifferentialFormValue{1,3}
@test isequal(Symbolics.simplify(dωb.data[1] - λ₂), Symbolics.simplify(Symbolics.Num(0)))
@test isequal(Symbolics.simplify(dωb.data[2] - λ₁), Symbolics.simplify(Symbolics.Num(0)))
@test isequal(Symbolics.simplify(dωb.data[3]), Symbolics.simplify(Symbolics.Num(0)))

ω1b = DifferentialFormValue{1,3}((λ₁*λ₂, λ₁*λ₃, λ₂*λ₃))
ddωb = exterior_derivative(exterior_derivative(ω1b, λ), λ)
@test all(isequal(Symbolics.simplify(c), Symbolics.simplify(Symbolics.Num(0)))
          for c in ddωb.data)

# ── Symbolic Koszul operator ─────────────────────────────────────────────────

# κ on a 1-form
ω1 = DifferentialFormValue{1,2}((x¹, x²))    # x¹ dx¹ + x² dx²
κω1 = koszul(v2, ω1)
@test typeof(κω1) <: DifferentialFormValue{0,2}
@test isequal(Symbolics.simplify(κω1.data[1]), Symbolics.simplify(x¹^2 + x²^2))

# κ on a 2-form
ω2 = DifferentialFormValue{2,2}((x¹,))       # x¹ dx¹∧dx²
κω2 = koszul(v2, ω2)
@test typeof(κω2) <: DifferentialFormValue{1,2}
@test isequal(Symbolics.simplify(κω2.data[1]), Symbolics.simplify(-x¹*x²))
@test isequal(Symbolics.simplify(κω2.data[2]), Symbolics.simplify(x¹^2))

# ── Homotopy identity (dκ + κd)(ω) = (r+K)ω  on P_r Λ^K ─────────────────────

function check_homotopy(ω::DifferentialFormValue{K,D}, r, x) where {K,D}
  dω  = exterior_derivative(ω)
  κdω = K < D ? koszul(x, dω) : zero(ω)   # κd (κ on (K+1)-form, only if K < D)
  κω  = K >= 1 ? koszul(x, ω) : ω         # κω (only if K ≥ 1, else 0 by convention)
  dκω = K >= 1 ? exterior_derivative(κω) : zero(ω)
  lhs = dκω + κdω
  rhs = (r + K) * ω
  ntuple(i -> isequal(Symbolics.simplify(lhs.data[i] - rhs.data[i]),
                      Symbolics.simplify(zero(x¹))), length(ω.data))
end

# K=0, D=2: only κd term (κ undefined on 0-forms, convention: 0)
# ω = x¹² (homogeneous deg 2) — κd(ω) should equal 2ω
ω_0 = DifferentialFormValue{0,2}((x¹^2,))
κdω_0 = koszul(v2, exterior_derivative(ω_0))
@test isequal(Symbolics.simplify(κdω_0.data[1] - 2*x¹^2), Symbolics.simplify(zero(x¹)))

# K=1, D=2, r=1: ω = x¹ dx¹ + x² dx²
@test all(check_homotopy(DifferentialFormValue{1,2}((x¹, x²)), 1, v2))

# K=1, D=2, r=1: ω = x¹ dx² (off-diagonal)
@test all(check_homotopy(DifferentialFormValue{1,2}((Symbolics.Num(0), x¹)), 1, v2))

# K=1, D=2, r=0: ω = dx¹ (constant 1-form)
@test all(check_homotopy(DifferentialFormValue{1,2}((Symbolics.Num(1), Symbolics.Num(0))), 0, v2))

# ── Symbolic codifferential ──────────────────────────────────────────────────

# In 2D on 1-forms with ⋆dx¹ = dx², ⋆dx² = -dx¹ and sign (-1)^{D(K-1)+1} = -1:
# δω = -(∂ω₁/∂x¹ + ∂ω₂/∂x²).  For ω₁=x¹²x², ω₂=x¹x²²: δω = -4x¹x²
f1 = DifferentialFormValue{1,2}((x¹^2 * x², x¹ * x²^2))
δf1 = codifferential(f1)
@test typeof(δf1) <: DifferentialFormValue{0,2}
@test isequal(Symbolics.simplify(δf1.data[1]), Symbolics.simplify(-4*x¹*x²))

# ── Lie derivative (Cartan formula) ──────────────────────────────────────────

# On a 0-form: L_v f = ι_v df = v¹ ∂f/∂x¹ + v² ∂f/∂x²
f0   = DifferentialFormValue{0,2}((x¹^2 + x²,))
v_e1 = VectorValue(Symbolics.Num(1), Symbolics.Num(0))
Lv_f0 = lie_derivative(v_e1, f0)
@test typeof(Lv_f0) <: DifferentialFormValue{0,2}
@test isequal(Symbolics.simplify(Lv_f0.data[1]), Symbolics.simplify(2*x¹))

# On a 1-form: ω = x¹dx¹ + x²dx², v = e₁
# ι_{e₁} ω = x¹, d(x¹) = dx¹; dω = 0  →  L_{e₁} ω = dx¹ = (1, 0)
ω_1 = DifferentialFormValue{1,2}((x¹, x²))
Lv_ω1 = lie_derivative(v_e1, ω_1)
@test typeof(Lv_ω1) <: DifferentialFormValue{1,2}
@test isequal(Symbolics.simplify(Lv_ω1.data[1]), Symbolics.simplify(Symbolics.Num(1)))
@test isequal(Symbolics.simplify(Lv_ω1.data[2]), Symbolics.simplify(Symbolics.Num(0)))

# ── print_forms: symbolic display of the PΛ bases ────────────────────────────

for (make, header) in ((RotatingPΛBasis, "RotatingPΛBasis"),
                       (TrimmedPΛBasis,  "TrimmedPΛBasis"),
                       ((V,T,r) -> BarycentricPΛBasis(V,T,r,1), "BarycentricPΛBasis"),
                       ((V,T,r) -> BarycentricPΛBasis(V,T,r,1; flavor=:BMM), "BarycentricPΛBasis"),
                       ((V,T,r) -> BarycentricPΛBasis(V,T,r,2), "BarycentricPΛBasis"),
                       ((V,T,r) -> BarycentricPmΛBasis(V,T,r,1), "BarycentricPmΛBasis"),
                       ((V,T,r) -> BarycentricPmΛBasis(V,T,r,1; flavor=:BMM), "BarycentricPmΛBasis"),
                       ((V,T,r) -> BarycentricPmΛBasis(V,T,r,2), "BarycentricPmΛBasis"))
  b   = make(Val(2), Float64, 2)
  buf = IOBuffer()
  print_forms(b, buf)
  out = String(take!(buf))
  @test occursin(header, out)
  @test occursin("dim = $(length(b))", out)
  # forms are printed on the ambient barycentric coframe, not the Cartesian one
  @test occursin("dλ", out)
  # one line per basis function
  @test count(==('['), out) >= length(b)
end

# The :BMM direction forms are the ψ of RotatingPΛBasis, so the two bases print
# the same ambient forms in the same order. r must be ≥ 3 for ψ to differ from
# the AFW direction form at all.
function form_lines(b)
  buf = IOBuffer()
  print_forms(b, buf)
  [ m.captures[1] for m in eachmatch(r"α=\(.*?\)\s+(.*)", String(take!(buf))) ]
end

let r = 3
  bmm = form_lines(BarycentricPΛBasis(Val(2), Float64, r, 1; flavor=:BMM))
  afw = form_lines(BarycentricPΛBasis(Val(2), Float64, r, 1))
  rot = form_lines(RotatingPΛBasis(Val(2), Float64, r))
  @test length(bmm) == length(rot) > 0
  @test bmm == rot
  @test afw != rot
end

# Likewise the :BMM P⁻ forms are the ϕ of TrimmedPΛBasis, scaled by the same
# bare monomial. Here the two bases enumerate their bubbles in a different order
# from r ≥ 2, so the printed lines match as sets rather than in sequence.
let r = 3
  bmm = form_lines(BarycentricPmΛBasis(Val(2), Float64, r, 1; flavor=:BMM))
  afw = form_lines(BarycentricPmΛBasis(Val(2), Float64, r, 1))
  trm = form_lines(TrimmedPΛBasis(Val(2), Float64, r))
  @test length(bmm) == length(trm) > 0
  @test sort(bmm) == sort(trm)
  @test sort(afw) != sort(trm)
end

end # module
