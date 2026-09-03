module DifferentialFormValuesTests
# Construction and pointwise algebra of DifferentialFormValue:
# zero/+/-/*, wedge product, interior product ι, flat Euclidean Hodge star.

using Gridap.TensorValues
using Test

# ── Construction for all (K,D) up to D=4 ─────────────────────────────────────

DifferentialFormValue{0,2}((1,))
DifferentialFormValue{1,2}((1,2))
DifferentialFormValue{2,2}((1,))

DifferentialFormValue{0,3}((1,))
DifferentialFormValue{1,3}((1,2,3))
DifferentialFormValue{2,3}((1,2,3))
DifferentialFormValue{3,3}((1,))

DifferentialFormValue{0,4}((1,))
DifferentialFormValue{1,4}((1,2,3,4))
DifferentialFormValue{2,4}((1,2,3,4,5,6))
DifferentialFormValue{3,4}((1,2,3,4))
DifferentialFormValue{4,4}((1,))

# ── zero, +, -, * ────────────────────────────────────────────────────────────

a3 = DifferentialFormValue{1,3}((1.0, 2.0, 3.0))
b3 = DifferentialFormValue{1,3}((4.0, 5.0, 6.0))

@test zero(a3).data == (0.0, 0.0, 0.0)
@test (a3 + b3).data == (5.0, 7.0, 9.0)
@test (b3 - a3).data == (3.0, 3.0, 3.0)
@test (-a3).data   == (-1.0, -2.0, -3.0)
@test (2.0 * a3).data == (2.0, 4.0, 6.0)
@test (a3 * 3.0).data == (3.0, 6.0, 9.0)

# s*(a+b) == s*a + s*b
s = 5.0
@test (s * (a3 + b3)).data == (s*a3 + s*b3).data

# ── Wedge product ────────────────────────────────────────────────────────────

# ∧ of 0-form and 1-forms in 2D
D = 2
a = DifferentialFormValue{0,D}((4,))
b = DifferentialFormValue{1,D}((1,1))
c = a ∧ b
@test typeof(c) <: DifferentialFormValue{1,D}
@test c.data == (4.0,4.0)

# ∧ of 1-forms in 2D
a = DifferentialFormValue{1,D}((0,1))
b = DifferentialFormValue{1,D}((1,0))
c = a ∧ b
@test typeof(c) <: DifferentialFormValue{2,D}
@test c.data == (-1.0,)

a = DifferentialFormValue{1,D}((1,0))
b = DifferentialFormValue{1,D}((0,1))
@test (a ∧ b).data == (1.0,)

a = DifferentialFormValue{1,D}((1,0))
b = DifferentialFormValue{1,D}((1,0))
@test (a ∧ b).data == (0.0,)

a = DifferentialFormValue{1,D}((0,1))
b = DifferentialFormValue{1,D}((0,1))
@test (a ∧ b).data == (0.0,)

a = DifferentialFormValue{1,D}((1,1))
b = DifferentialFormValue{1,D}((1,1))
@test (a ∧ b).data == (0.0,)

a = DifferentialFormValue{1,D}((4,1))
b = DifferentialFormValue{1,D}((1,1))
@test (a ∧ b).data == (3.0,)

# ∧ of 1-form and 2-form in 2D: K1+K2 = 3 > D = 2, result is the zero 3-form
a = DifferentialFormValue{1,D}((4,1))
b = DifferentialFormValue{2,D}((1,))
c = a ∧ b
@test typeof(c) <: DifferentialFormValue{3,D}
@test c.data == ()

# ∧ of 0-form and 2-form in 2D
a = DifferentialFormValue{0,D}((4.0,))
b = DifferentialFormValue{2,D}((1.0,))
c = a ∧ b
@test typeof(c) <: DifferentialFormValue{2,D}
@test c.data == (4.0,)

# 3D wedges
D = 3
a = DifferentialFormValue{1,D}((0,0,1))
b = DifferentialFormValue{1,D}((1,0,0))
c = a ∧ b
@test typeof(c) <: DifferentialFormValue{2,D}
@test c.data == (0.0,-1.0,0.0)

a = DifferentialFormValue{1,D}((1,0,0))
b = DifferentialFormValue{1,D}((1,1,0))
@test (a ∧ b).data == (1.0,0.0,0.0)

a = DifferentialFormValue{1,D}((1,0,0))
b = DifferentialFormValue{2,D}((1,1,0))
c = a ∧ b
@test typeof(c) <: DifferentialFormValue{3,D}
@test c.data == (0.0,)

a = DifferentialFormValue{1,D}((0,0,1))
b = DifferentialFormValue{2,D}((1,0,0))
@test (a ∧ b).data == (1.0,)

a = DifferentialFormValue{0,D}((4.0,))
b = DifferentialFormValue{2,D}((1.0,0.0,0.0))
c = a ∧ b
@test typeof(c) <: DifferentialFormValue{2,D}
@test c.data == (4.0,0.0,0.0)

# ── Interior product ι ───────────────────────────────────────────────────────

D2 = 2
e1 = VectorValue(1.0, 0.0)
e2 = VectorValue(0.0, 1.0)

dx1_2d = DifferentialFormValue{1,D2}((1.0, 0.0))
dx2_2d = DifferentialFormValue{1,D2}((0.0, 1.0))

# ι_eᵢ(dxʲ) = δᵢʲ  (0-form = scalar)
r11 = ι(e1, dx1_2d);  @test typeof(r11) <: DifferentialFormValue{0,D2};  @test r11.data == (1.0,)
r22 = ι(e2, dx2_2d);  @test r22.data == (1.0,)
r12 = ι(e1, dx2_2d);  @test r12.data == (0.0,)
r21 = ι(e2, dx1_2d);  @test r21.data == (0.0,)

# ι_e₁(dx¹∧dx²) = dx², ι_e₂(dx¹∧dx²) = -dx¹
vol2 = DifferentialFormValue{2,D2}((1.0,))
rv1  = ι(e1, vol2);  @test typeof(rv1) <: DifferentialFormValue{1,D2};  @test rv1.data == (0.0, 1.0)
rv2  = ι(e2, vol2);  @test rv2.data == (-1.0, 0.0)

# In 3D: ι_e₁(dx¹∧dx²∧dx³) = dx²∧dx³
# 2-combinations in 3D ordered: (1,2),(1,3),(2,3)
D3   = 3
e13  = VectorValue(1.0, 0.0, 0.0)
vol3 = DifferentialFormValue{3,D3}((1.0,))
r3   = ι(e13, vol3)
@test typeof(r3) <: DifferentialFormValue{2,D3}
@test r3.data == (0.0, 0.0, 1.0)   # dx²∧dx³ is the third 2-combination

# ── Hodge star (flat Euclidean) ──────────────────────────────────────────────

# 3D: ⋆dx¹ = dx²∧dx³,  ⋆dx² = -dx¹∧dx³,  ⋆dx³ = dx¹∧dx²
dx1_3d = DifferentialFormValue{1,D3}((1.0, 0.0, 0.0))
dx2_3d = DifferentialFormValue{1,D3}((0.0, 1.0, 0.0))
dx3_3d = DifferentialFormValue{1,D3}((0.0, 0.0, 1.0))

@test hodge_star(dx1_3d).data == (0.0,  0.0,  1.0)    # dx²∧dx³
@test hodge_star(dx2_3d).data == (0.0, -1.0,  0.0)    # -dx¹∧dx³
@test hodge_star(dx3_3d).data == (1.0,  0.0,  0.0)    # dx¹∧dx²

# ⋆∘⋆ = +1 for K=1, D=3 ((-1)^{K(D-K)} = (-1)^2 = +1)
for ω in [dx1_3d, dx2_3d, dx3_3d]
  @test hodge_star(hodge_star(ω)).data == ω.data
end

# 0-form → 3-form
f0 = DifferentialFormValue{0,D3}((5.0,))
sf0 = hodge_star(f0)
@test typeof(sf0) <: DifferentialFormValue{3,D3}
@test sf0.data == (5.0,)

# 3-form → 0-form
svol3 = hodge_star(vol3)
@test typeof(svol3) <: DifferentialFormValue{0,D3}
@test svol3.data == (1.0,)

# 2D: ⋆ maps 1-forms to 1-forms (Kc = D-K = 2-1 = 1)
# ⋆dx¹ = dx², ⋆dx² = -dx¹  (standard 2D orientation)
sdx1_2d = hodge_star(dx1_2d)
sdx2_2d = hodge_star(dx2_2d)
@test typeof(sdx1_2d) <: DifferentialFormValue{1,D2}
@test sdx1_2d.data == (0.0, 1.0)    # 0*dx¹ + 1*dx²
@test sdx2_2d.data == (-1.0, 0.0)   # -1*dx¹ + 0*dx²

# ⋆∘⋆ = -1 for K=1, D=2 ((-1)^{K(D-K)} = (-1)^1 = -1)
for ω in [dx1_2d, dx2_2d]
  @test hodge_star(hodge_star(ω)).data == (-1.0 * ω).data
end

# ── Display: coframe labelling ───────────────────────────────────────────────

_str(io_props, ω) = sprint((io, x) -> show(io, MIME("text/plain"), x), ω;
                           context = io_props)

ω_show = DifferentialFormValue{1,D3}((1.0, 2.0, 3.0))
@test _str(:coordinates => :cartesian,   ω_show) == "1.0 dx¹ + 2.0 dx² + 3.0 dx³"
@test _str(:coordinates => :barycentric, ω_show) == "1.0 dλ¹ + 2.0 dλ² + 3.0 dλ³"
# Cartesian labelling is the default when the property is absent
@test sprint((io, x) -> show(io, MIME("text/plain"), x), ω_show) ==
      _str(:coordinates => :cartesian, ω_show)
# An unrecognised coframe is an error, not a silent fallback
@test_throws ArgumentError _str(:coordinates => :something_else, ω_show)
@test_throws ArgumentError _str(:coordinates => true, ω_show)

# Wedge basis elements are labelled on both coframes
ω2_show = DifferentialFormValue{2,D3}((1.0, 0.0, 0.0))
@test occursin("dx¹ ∧ dx²", _str(:coordinates => :cartesian,   ω2_show))
@test occursin("dλ¹ ∧ dλ²", _str(:coordinates => :barycentric, ω2_show))

end # module
