module DifferentialFormValuesTests
# Construction and pointwise algebra of DifferentialFormValue:
# zero/+/-/*, wedge product, interior product ι, flat Euclidean Hodge star,
# inner product of forms.

using Gridap.TensorValues
using StaticArrays
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

# ── Constructor family ───────────────────────────────────────────────────────

# Homogeneous NTuple, with and without an explicit component type
@test isa(DifferentialFormValue{1,3}((1,2,3)),            DifferentialFormValue{1,3,Int,3})
@test isa(DifferentialFormValue{1,3,Float64}((1,2,3)),    DifferentialFormValue{1,3,Float64,3})
@test isa(DifferentialFormValue{1,3,Float64,3}((1,2,3)),  DifferentialFormValue{1,3,Float64,3})
@test isa(DifferentialFormValue{1,3,Float64,3}((1.,2.,3.)), DifferentialFormValue{1,3,Float64,3})

# Heterogeneous Tuple: promoted, or converted to the requested type
@test isa(DifferentialFormValue{1,3}((1,2.0,3)),          DifferentialFormValue{1,3,Float64,3})
@test isa(DifferentialFormValue{1,3,Int}((1,2.0,3)),      DifferentialFormValue{1,3,Int,3})
@test isa(DifferentialFormValue{1,3,Int,3}((1,2.0,3)),    DifferentialFormValue{1,3,Int,3})

# Vararg
@test DifferentialFormValue{1,3}(1,2.0,3)       == DifferentialFormValue{1,3}((1.0,2.0,3.0))
@test DifferentialFormValue{1,3,Int}(1,2.0,3)   == DifferentialFormValue{1,3}((1,2,3))
@test DifferentialFormValue{1,3,Int,3}(1,2.0,3) == DifferentialFormValue{1,3}((1,2,3))

# A single scalar is a full 0-form and a full D-form
@test isa(DifferentialFormValue{0,3}(1),   DifferentialFormValue{0,3,Int,1})
@test isa(DifferentialFormValue{0,0}(1),   DifferentialFormValue{0,0,Int,1})
@test isa(DifferentialFormValue{3,3}(1.0), DifferentialFormValue{3,3,Float64,1})

# K > D: the space is trivial, so the value has no component
@test isa(DifferentialFormValue{4,3}(),          DifferentialFormValue{4,3,Int,0})
@test isa(DifferentialFormValue{4,3}(()),        DifferentialFormValue{4,3,Int,0})
@test isa(DifferentialFormValue{4,3,Float64}(),  DifferentialFormValue{4,3,Float64,0})
@test isa(DifferentialFormValue{4,3,Float64}(()),DifferentialFormValue{4,3,Float64,0})
@test isa(DifferentialFormValue{4,3,Int,0}(),    DifferentialFormValue{4,3,Int,0})
@test isa(DifferentialFormValue{4,3,Int,0}(()),  DifferentialFormValue{4,3,Int,0})

@test_throws AssertionError DifferentialFormValue{1,3}((1,2))
@test_throws AssertionError DifferentialFormValue{1,3,Int}((1,2))
@test_throws AssertionError DifferentialFormValue{1,3,Int,2}((1,2))
@test_throws AssertionError DifferentialFormValue{1,3,Int,2}(1,2)
@test_throws AssertionError DifferentialFormValue{0,3}()
@test_throws AssertionError DifferentialFormValue{0,3,Int}()

# ── Independent component count ──────────────────────────────────────────────

@test num_indep_components(DifferentialFormValue{0,2}) == 1
@test num_indep_components(DifferentialFormValue{2,2}) == 1
@test num_indep_components(DifferentialFormValue{1,3}) == 3
@test num_indep_components(DifferentialFormValue{2,3}) == 3
@test num_indep_components(DifferentialFormValue{4,3}) == 0

@test zero(DifferentialFormValue{2,2,Float64}) == DifferentialFormValue{2,2}((0.0,))
@test zero(DifferentialFormValue{2,3,Float64}) == DifferentialFormValue{2,3}((0.0,0.0,0.0))

# ── Conversion to a static array ─────────────────────────────────────────────

# The stored components expand into the full D^K antisymmetric tensor, with the
# sorted multi-indices carrying the components themselves
A = convert(SArray{Tuple{3,3},Int}, DifferentialFormValue{2,3}((1,2,3)))
@test A == [0 1 2; -1 0 3; -2 -3 0]
@test A == -transpose(A)
@test isa(convert(MArray{Tuple{3,3},Int}, DifferentialFormValue{2,3}((1,2,3))), MMatrix{3,3,Int,9})

@test eltype(convert(SArray{Tuple{3,3},Float64}, DifferentialFormValue{2,3}((1,2,3)))) == Float64

# A 1-form gives back its components, a 0-form a zero-dimensional array
@test convert(SArray{Tuple{3},Int}, DifferentialFormValue{1,3}((1,2,3))) == [1,2,3]
@test convert(SArray{Tuple{},Int}, DifferentialFormValue{0,3}((7,)))[] == 7

# The volume form is the Levi-Civita symbol scaled by its single component
A3 = convert(SArray{Tuple{3,3,3},Int}, DifferentialFormValue{3,3}((2,)))
@test (A3[1,2,3], A3[2,1,3], A3[3,1,2], A3[1,1,3]) == (2, -2, 2, 0)

# K > D: no stored component, but the tensor still has D^K (vanishing) entries
A4 = convert(SArray{Tuple{3,3,3,3},Int}, DifferentialFormValue{4,3}())
@test size(A4) == (3,3,3,3)
@test all(iszero, A4)

# The wedge of two 1-forms expands to the antisymmetric part of their outer product
a1 = DifferentialFormValue{1,4}((1,2,3,4))
b1 = DifferentialFormValue{1,4}((5,6,7,8))
va = convert(SArray{Tuple{4},Int}, a1)
vb = convert(SArray{Tuple{4},Int}, b1)
@test convert(SArray{Tuple{4,4},Int}, a1 ∧ b1) == va*transpose(vb) - vb*transpose(va)

# A shape that is not the form's own must not convert
@test_throws DimensionMismatch convert(SArray{Tuple{2,2},Int}, DifferentialFormValue{2,3}((1,2,3)))

# Internal conversion changes the component type, and is a no-op on its own type
@test convert(DifferentialFormValue{2,3,Float64}, DifferentialFormValue{2,3}((1,2,3))) ==
      DifferentialFormValue{2,3}((1.0,2.0,3.0))
ω23 = DifferentialFormValue{2,3}((1.0,2.0,3.0))
@test convert(DifferentialFormValue{2,3,Float64}, ω23) === ω23

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

# ── Hodge star signs against the full Levi-Civita matrix ─────────────────────

# The dense sign matrix ε_{I_n J_m}, over every pair of a K-combination I_n and
# a (D−K)-combination J_m of 1:D. `_hodge_star_signs` is the claim that this is
# zero off the anti-diagonal, where it holds ε_{I_n I_nᶜ}.
function reference_hodge_star_matrix(K::Int, D::Int)
  c_in  = sorted_combinations(D, K)
  c_out = sorted_combinations(D, D-K)
  M = zeros(Int, length(c_out), length(c_in))
  for (m, J) in enumerate(c_out), (n, I) in enumerate(c_in)
    perm = [I..., J...]
    length(unique(perm)) == D || continue
    M[m, n] = sorting_sign(perm...)
  end
  M
end

for D in 1:6, K in 0:D
  L = binomial(D, K)
  M = reference_hodge_star_matrix(K, D)
  s = TensorValues._hodge_star_signs(K, D)

  @test length(s) == L
  @test size(M) == (L, L)
  for m in 1:L, n in 1:L
    @test M[m,n] == (n == L+1-m ? s[n] : 0)
  end
end

# ── Inner product of forms ───────────────────────────────────────────────────

# The order-K tensor holding the expanded components of a form
_full(ω::DifferentialFormValue) = MultiValue(SArray(ω))

# An inverse metric with no zero minor, so every term of the raising sum counts
g_inv_fi = SymTensorValue{D3}(2.0, -0.5, 0.3, 1.5, 0.25, 3.0)
sdg_fi   = 1 / sqrt(det(g_inv_fi))

for (α, β) in [(DifferentialFormValue{0,D3}((7.0,)),
                DifferentialFormValue{0,D3}((3.0,))),
               (DifferentialFormValue{1,D3}((1.0, 2.0, 3.0)),
                DifferentialFormValue{1,D3}((4.0, 5.0, 6.0))),
               (DifferentialFormValue{2,D3}((1.0, 2.0, 3.0)),
                DifferentialFormValue{2,D3}((4.0, 5.0, 6.0))),
               (DifferentialFormValue{3,D3}((2.0,)),
                DifferentialFormValue{3,D3}((5.0,)))]
  K = length(size(α))

  @test α ⨟ β == form_inner(α, β)
  @test form_inner(α, β) == sum(α.data .* β.data)
  @test form_inner(α, β) ≈ form_inner(α, β, one(SymTensorValue{D3,Float64}))
  # α ∧ ⋆_g β = √det(g) (α|β) dx¹∧…∧dx^D, the property that defines (α|β)
  @test vol_coeff(α ∧ hodge_star(β, g_inv_fi, sdg_fi)) ≈ sdg_fi * form_inner(α, β, g_inv_fi)
  # the full contraction of the expanded components drops the 1/K! of (α|β)
  @test contracted_product(Val(K), α, β) ≈ factorial(K) * form_inner(α, β)
  # summing the binomial(D,K) independent components gives the same as
  # contracting all D^K of them
  @test inner(α, β) ≈ contracted_product(Val(K), α, β)
  if K > 0 # MultiValue has no constructor from the 0-dimensional array of a 0-form
    @test contracted_product(Val(K), _full(α), _full(β)) ≈ factorial(K) * form_inner(α, β)
  end
end

# A 0-form has no index to contract, it only enters the tensor product
ω0_fi = DifferentialFormValue{0,D3}((2.0,))
ω1_fi = DifferentialFormValue{1,D3}((1.0, 2.0, 3.0))
@test outer(ω0_fi, ω1_fi) === VectorValue{D3}(2.0, 4.0, 6.0)
@test_throws DimensionMismatch contracted_product(Val(1), ω0_fi, ω1_fi)

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

# ── Scalar types: mixed arguments and degenerate shapes ──────────────────────

# The eltype comes from the arithmetic on the components, so it follows the
# usual promotion of the scalar types of the arguments.
a_i  = DifferentialFormValue{1,D3,Int}((1,2,3))
b_f  = DifferentialFormValue{1,D3,Float32}((1f0,2f0,3f0))
v_i  = VectorValue{D3,Int}(1,2,3)
g_f  = SymTensorValue{D3,Float32}(1f0,0f0,0f0,1f0,0f0,1f0)
Jt_f = TensorValue{D3,2,Float32}(1f0,0f0,0f0,0f0,1f0,0f0)

@test eltype(a_i ∧ b_f)                 === Float32
@test eltype(ι(v_i, b_f))               === Float32
@test eltype(hodge_star(a_i))           === Int
@test eltype(hodge_star(a_i, g_f, 1.0)) === Float64
@test eltype(flat(v_i, g_f))            === Float32
@test eltype(sharp(b_f, SymTensorValue{D3,Int}(1,0,0,1,0,1))) === Float32
@test eltype(pullback(a_i, Jt_f))       === Float32
@test eltype(outer(v_i, b_f))           === Float32
@test typeof(form_inner(a_i, b_f))      === Float32
@test typeof(form_inner(a_i, b_f, g_f)) === Float32
@test typeof(apply_form(a_i, v_i))      === Int
@test eltype(grad_to_2form(TensorValue{D3,D3,Int}(1,2,3,4,5,6,7,8,9))) === Int

# a K-form takes exactly K arguments, each a size (D,) tensor
@test_throws ErrorException apply_form(a_i, v_i, v_i)
@test_throws ErrorException apply_form(a_i)
@test_throws ErrorException apply_form(a_i, (1,2,3))
@test_throws ErrorException apply_form(a_i, [1,2,3])
@test_throws ErrorException apply_form(a_i, VectorValue{2,Int}(1,2))
@test_throws ErrorException apply_form(a_i, TensorValue{1,D3,Int}(1,2,3))

# Λᴷ is trivial for K > D. Operations on such a form must keep the promoted
# scalar type rather than fall back to the eltype of the empty tuple.
triv = DifferentialFormValue{4,D3,Float32}(())

# ι lands in the nontrivial Λ³: a zero form, not an empty one
@test ι(v_i, triv) === zero(DifferentialFormValue{3,D3,Float32})
@test a_i ∧ triv   === zero(DifferentialFormValue{5,D3,Float32})
@test pullback(triv, Jt_f) === zero(DifferentialFormValue{4,2,Float32})
@test apply_form(triv, v_i, v_i, v_i, v_i) === zero(Float32)
# the vectors take part in the promotion even though no term survives
@test (@inferred apply_form(triv, v_i, v_i, v_i, VectorValue{D3,Float64}(1,2,3))) === zero(Float64)
@test apply_form(DifferentialFormValue{4,D3,Int}(()), v_i, v_i, v_i, v_i) === zero(Int)
@test apply_form(DifferentialFormValue{1,0,Float32}(()), VectorValue{0,Float32}()) === zero(Float32)
@test form_inner(triv, DifferentialFormValue{4,D3,Int}(())) === zero(Float32)
@test form_inner(triv, DifferentialFormValue{4,D3,Int}(()), g_f) === zero(Float32)
# the conversions name the eltype, so it survives an empty container
@test to_Kform(VectorValue{0,Float32}(()), Val(4), Val(D3)) ===
      zero(DifferentialFormValue{4,D3,Float32})
@test from_Kform(triv) === zero(VectorValue{0,Float32})
@test to_1form(VectorValue{0,Float32}()) === zero(DifferentialFormValue{1,0,Float32})
@test from_1form(DifferentialFormValue{1,0,Float32}(())) === zero(VectorValue{0,Float32})

# ⋆ of such a form would be a Λ^(D-K) with D-K < 0, which is not a type
@test_throws ErrorException hodge_star(triv)
@test_throws ErrorException hodge_star(triv, SymTensorValue{D3,Int}(1,0,0,1,0,1), 1.0)

end # module
