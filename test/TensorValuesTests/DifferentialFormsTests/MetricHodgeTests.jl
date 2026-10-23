module MetricHodgeTests
# Metric-aware Hodge star and musical operators (flat/sharp).
# Uses the CubedSphere metric (gnomonic projection of a cube face onto the
# sphere) as a concrete non-trivial metric, plus a diagonal metric with
# analytically exact values.

using Gridap.TensorValues
using LinearAlgebra
using Test

# ── CubedSphere metric helpers ────────────────────────────────────────────────
#
# Gnomonic projection from a cube face to the sphere of radius r.
# With a = tan α, b = tan β, sa = 1+a², sb = 1+b², ρ² = 1+a²+b²:
#
#   g₁₁ = r²·sa²·sb / ρ⁴,   g₁₂ = −r²·a·b·sa·sb / ρ⁴,   g₂₂ = r²·sa·sb² / ρ⁴
#   det(g) = r⁴·sa²·sb² / ρ⁶   →   √det(g) = r²·sa·sb / ρ³

function _csphere_metric(r, x)
  a, b   = tan(x[1]), tan(x[2])
  sa, sb = 1 + a^2, 1 + b^2
  ρ4     = (1 + a^2 + b^2)^2
  c      = r^2 * sa * sb / ρ4
  SymTensorValue{cov2,Float64}(c*sa, -c*a*b, c*sb)
end

function _csphere_inv_metric(r, x)
  a, b   = tan(x[1]), tan(x[2])
  sa, sb = 1 + a^2, 1 + b^2
  ρ2     = 1 + a^2 + b^2
  c      = ρ2 / (r^2 * sa * sb)
  SymTensorValue{con2,Float64}(c*sb, c*a*b, c*sa)
end

function _csphere_sqrt_det_g(r, x)
  a, b   = tan(x[1]), tan(x[2])
  sa, sb = 1 + a^2, 1 + b^2
  ρ3     = (1 + a^2 + b^2)^(3/2)
  r^2 * sa * sb / ρ3
end

# ── Verify metric self-consistency ───────────────────────────────────────────

r  = 1.0
x0 = VectorValue(π/6, π/4)   # (30°, 45°) — generic asymmetric point
g     = _csphere_metric(r, x0)
g_inv = _csphere_inv_metric(r, x0)
sdg   = _csphere_sqrt_det_g(r, x0)

# g · g⁻¹ ≈ I₂
I2 = TensorValue{2,2,Float64,4}(1.0, 0.0, 0.0, 1.0)
g_full     = TensorValue{2,2,Float64,4}(g[1,1],     g[2,1],     g[1,2],     g[2,2])
g_inv_full = TensorValue{2,2,Float64,4}(g_inv[1,1], g_inv[2,1], g_inv[1,2], g_inv[2,2])
@test norm(g_full ⋅ g_inv_full - I2) < 1e-12

@test abs(g[1,1]*g[2,2] - g[1,2]^2 - sdg^2) < 1e-12

# ── Part 1: Diagonal metric (analytical exact values) ────────────────────────
#
# g = diag(p, q),  g_inv = diag(1/p, 1/q),  √det g = √(pq)
# ⋆_g dx¹ = √(p/q) dx²    ⋆_g dx² = −√(q/p) dx¹

p, q       = 2.0, 3.0
g_diag     = SymTensorValue{cov2,Float64,3}(p, 0.0, q)
g_inv_diag = inv(g_diag) # con2
sdg_diag   = sqrt(p*q)

dx1 = CoVectorValue{2}((1.0, 0.0))
dx2 = CoVectorValue{2}((0.0, 1.0))

star_dx1 = hodge_star(dx1, g_inv_diag, sdg_diag)
star_dx2 = hodge_star(dx2, g_inv_diag, sdg_diag)

@test typeof(star_dx1) <: CoVectorValue{cov2}
@test collect(star_dx1.data) ≈ [0.0, sqrt(q/p)]    atol=1e-12   # √(q/p) dx²
@test collect(star_dx2.data) ≈ [-sqrt(p/q), 0.0]   atol=1e-12   # -√(p/q) dx¹

# ⋆_g ∘ ⋆_g = (-1)^{K(D-K)} = -1 for K=1, D=2
@test collect(hodge_star(star_dx1, g_inv_diag, sdg_diag).data) ≈ collect((-1.0 * dx1).data)  atol=1e-12
@test collect(hodge_star(star_dx2, g_inv_diag, sdg_diag).data) ≈ collect((-1.0 * dx2).data)  atol=1e-12

# 0-form: ⋆_g f = f √det(g) vol_coord
f0    = ExteriorFormValue{0,2}((5.0,))
sf0_d = hodge_star(f0, g_inv_diag, sdg_diag)
@test typeof(sf0_d) <: ExteriorFormValue{2,cov2}
@test collect(sf0_d.data) ≈ [5.0 * sdg_diag]   atol=1e-12

# 2-form: ⋆_g vol_coord = 1/√det(g)
vol_coord = ExteriorFormValue{2,2}((1.0,))
svol_d    = hodge_star(vol_coord, g_inv_diag, sdg_diag)
@test typeof(svol_d) <: ExteriorFormValue{0,cov2}
@test collect(svol_d.data) ≈ [1.0/sdg_diag]   atol=1e-12

# Supplying `g` instead of `g_inv`, and omitting `√det g`, are equivalent: the
# missing data is derived from whichever tensor is given.
for ω in [dx1, dx2, f0, vol_coord]
  ref = hodge_star(ω, g_inv_diag, sdg_diag)
  @test collect(hodge_star(ω, g_inv_diag).data)        ≈ collect(ref.data)  atol=1e-12
  @test collect(hodge_star(ω, g_diag, sdg_diag).data)  ≈ collect(ref.data)  atol=1e-12
  @test collect(hodge_star(ω, g_diag).data)            ≈ collect(ref.data)  atol=1e-12
end

# ── Part 2: CubedSphere metric at a generic point ─────────────────────────────

# ⋆_g ∘ ⋆_g = -1 for K=1, D=2
for ω in [dx1, dx2]
  sω  = hodge_star(ω,  g_inv, sdg)
  ssω = hodge_star(sω, g_inv, sdg)
  @test collect(ssω.data) ≈ collect((-1.0 * ω).data)   atol=1e-10
end

# Generic 1-form
omega_gen = CoVectorValue{2}((0.7, -0.3))
s1 = hodge_star(omega_gen, g_inv, sdg)
s2 = hodge_star(s1,        g_inv, sdg)
@test collect(s2.data) ≈ collect((-1.0 * omega_gen).data)   atol=1e-10

# 0-form ⋆_g: result is a 2-form with coefficient f·√det(g)
f0_cs  = ExteriorFormValue{0,2}((3.0,))
sf0_cs = hodge_star(f0_cs, g_inv, sdg)
@test typeof(sf0_cs) <: ExteriorFormValue{2,cov2}
@test collect(sf0_cs.data) ≈ [3.0 * sdg]   atol=1e-10

# Same three routes on a non-diagonal metric
for ω in [dx1, dx2, omega_gen, f0_cs]
  ref = hodge_star(ω, g_inv, sdg)
  @test collect(hodge_star(ω, g_inv).data)   ≈ collect(ref.data)  atol=1e-10
  @test collect(hodge_star(ω, g, sdg).data)  ≈ collect(ref.data)  atol=1e-10
  @test collect(hodge_star(ω, g).data)       ≈ collect(ref.data)  atol=1e-10
end

# ── Part 3: Musical operators ─────────────────────────────────────────────────

e1 = VectorValue(1.0, 0.0)
e2 = VectorValue(0.0, 1.0)

# Flat metric: flat/sharp are identities
g_flat     = one(SymTensorValue{cov2,Float64})
g_inv_flat = inv(g_flat)

@test flat(e1, g_flat).data == (1.0, 0.0)
@test flat(e2, g_flat).data == (0.0, 1.0)
@test sharp(dx1, g_inv_flat).data == (1.0, 0.0)
@test sharp(dx2, g_inv_flat).data == (0.0, 1.0)

# Diagonal metric: flat(eᵢ, g) lowers with g
@test collect(flat(e1, g_diag).data) ≈ [p, 0.0]
@test collect(flat(e2, g_diag).data) ≈ [0.0, q]

# Round-trip: sharp(flat(v, g), g_inv) = v
for v in [e1, e2, VectorValue(1.5, -0.7)]
  v_rt = sharp(flat(v, g_diag), g_inv_diag)
  @test collect(v_rt.data) ≈ collect(v.data)   atol=1e-12
end

# Round-trip: flat(sharp(ω, g_inv), g) = ω
for ω in [dx1, dx2, CoVectorValue{2}((0.4, -0.9))]
  ω_rt = flat(sharp(ω, g_inv_diag), g_diag)
  @test collect(ω_rt.data) ≈ collect(ω.data)   atol=1e-12
end

# CubedSphere: sharp ∘ flat = id
for v in [e1, e2, VectorValue(0.3, -0.8)]
  v_rt = sharp(flat(v, g), g_inv)
  @test collect(v_rt.data) ≈ collect(v.data)   atol=1e-10
end

# CubedSphere: flat ∘ sharp = id
for ω in [dx1, dx2, omega_gen]
  ω_rt = flat(sharp(ω, g_inv), g)
  @test collect(ω_rt.data) ≈ collect(ω.data)   atol=1e-10
end

# ── Part 4: L²-inner product via ω ∧ ⋆_g η ──────────────────────────────────
# For K=1, D=2 with diagonal metric g = diag(p,q):
#   ⟨dx¹, dx¹⟩_g = g^{11} = 1/p  →  dx¹ ∧ ⋆_g dx¹ = (1/p)·√(pq) dx¹∧dx²
#   ⟨dx², dx²⟩_g = g^{22} = 1/q  →  dx² ∧ ⋆_g dx² = (1/q)·√(pq) dx¹∧dx²
#   ⟨dx¹, dx²⟩_g = 0

inner11 = (dx1 ∧ hodge_star(dx1, g_inv_diag, sdg_diag)).data[1]
inner22 = (dx2 ∧ hodge_star(dx2, g_inv_diag, sdg_diag)).data[1]
inner12 = (dx1 ∧ hodge_star(dx2, g_inv_diag, sdg_diag)).data[1]

@test inner11 ≈ sdg_diag / p   atol=1e-12    # √(pq)/p = √(q/p)
@test inner22 ≈ sdg_diag / q   atol=1e-12    # √(pq)/q = √(p/q)
@test abs(inner12) < 1e-12

# ── Part 5: metric arguments ─────────────────────────────────────────────────
# Every metric operation takes the same trailing `metric_args...`: a
# `MetricData` is accepted wherever `g`/`g_inv`[, `√det g`] is.

md_lazy  = parse_metric_data(g_diag)
md_eager = compute_metric_data(g_diag)
v_gen    = ContraVectorValue{con2}((1.5, -0.7))
for md in (md_lazy, md_eager)
  @test collect(hodge_star(dx1, md).data)          ≈ collect(hodge_star(dx1, g_inv_diag, sdg_diag).data)  atol=1e-12
  @test collect(flat(v_gen, md).data)              ≈ collect(flat(v_gen, g_diag).data)                    atol=1e-12
  @test collect(sharp(omega_gen, md).data)         ≈ collect(sharp(omega_gen, g_inv_diag).data)           atol=1e-12
  @test alt_inner(dx1, omega_gen, md)              ≈ alt_inner(dx1, omega_gen, g_inv_diag)                 atol=1e-12
  @test alt_inner(v_gen, v_gen, md)                ≈ alt_inner(v_gen, v_gen, g_diag)                       atol=1e-12
  # a non-alternating tensor is starred through 𝔸, here dx¹∧dx² (data is column-major)
  @test collect(hodge_star(TensorValue{cov2,cov2}(0.,-1.,1.,0.), md).data) ≈ collect(hodge_star(ExteriorFormValue{2,cov2}((1.,)), g_diag).data)  atol=1e-12
end

# The metric must have the dimension of the tensor it acts on, whichever form
# it is given in
g3   = SymTensorValue{cov3}(1.,0.,0.,1.,0.,1.)
md3  = parse_metric_data(g3)
dx1_3 = CoVectorValue{cov3}((1.,0.,0.))
@test_throws DimensionMismatch hodge_star(dx1, g3)
@test_throws DimensionMismatch hodge_star(dx1, md3)
@test_throws DimensionMismatch flat(v_gen, g3)
@test_throws DimensionMismatch flat(e1, md3)
@test_throws DimensionMismatch sharp(omega_gen, inv(g3), 1.0)
@test_throws DimensionMismatch sharp(omega_gen, md3)
@test_throws DimensionMismatch alt_inner(dx1, dx1, g3)
@test_throws DimensionMismatch alt_inner(dx1, dx1, md3)
@test hodge_star(dx1_3, md3) isa ExteriorFormValue{2,cov3}

# Missing metric information is an error, not a silent Euclidean fallback
@test_throws ArgumentError hodge_star(dx1)
@test_throws ArgumentError flat(v_gen)
@test_throws ArgumentError sharp(omega_gen)
@test_throws ArgumentError alt_inner(dx1, dx1)

# The tensors and the scalar are told apart by type, so they parse in any order
# and only what is given is stored
for args in ((g_diag, g_inv_diag, sdg_diag), (sdg_diag, g_inv_diag, g_diag), (g_inv_diag, sdg_diag, g_diag))
  md = parse_metric_data(args...)
  @test getfield(md, :g)          === g_diag
  @test getfield(md, :g_inv)      === g_inv_diag
  @test getfield(md, :signed_sqrt_g) === sdg_diag
end
md_rev = parse_metric_data(sdg_diag, g_inv_diag)
@test getfield(md_rev, :g)          === nothing
@test getfield(md_rev, :g_inv)      === g_inv_diag
@test getfield(md_rev, :signed_sqrt_g) === sdg_diag
@test collect(hodge_star(dx1, sdg_diag, g_inv_diag).data) ≈ collect(hodge_star(dx1, g_inv_diag, sdg_diag).data)  atol=1e-12
@test compute_metric_data(sdg_diag, g_diag).signed_sqrt_g === sdg_diag
@test compute_metric_data(g_diag, sdg_diag).signed_sqrt_g === sdg_diag
# eager computation from either tensor alone stores all three data
for md in (compute_metric_data(g_diag), compute_metric_data(g_inv_diag), compute_metric_data(parse_metric_data(g_inv_diag)))
  @test getfield(md, :g)     isa SymTensorValue{cov2}
  @test getfield(md, :g_inv) isa SymTensorValue{con2}
  @test getfield(md, :signed_sqrt_g) isa Float64
  @test collect(md.g.data) ≈ collect(g_diag.data)  atol=1e-12
  @test collect(md.g_inv.data) ≈ collect(g_inv_diag.data)  atol=1e-12
  @test md.signed_sqrt_g ≈ sdg_diag  atol=1e-12
end

# Each datum at most once, tensors must carry a variance, anything else is rejected
@test_throws ArgumentError parse_metric_data(g_diag, g_diag)
@test_throws ArgumentError parse_metric_data(g_inv_diag, sdg_diag, g_inv_diag)
@test_throws ArgumentError parse_metric_data(g_diag, sdg_diag, 2.0)
@test_throws ArgumentError parse_metric_data(SymTensorValue{2}(1.,0.,1.))
@test_throws ArgumentError parse_metric_data(TensorValue{cov2,cov2}(1.,0.,0.,1.))
@test_throws ArgumentError parse_metric_data(g_diag, "signed_sqrt_g")
@test_throws ArgumentError parse_metric_data(md_lazy, sdg_diag)
@test_throws ArgumentError parse_metric_data(I, sdg_diag)
@test_throws ArgumentError hodge_star(dx1, md_lazy, sdg_diag)
@test_throws ArgumentError flat(v_gen, g_diag, TensorValue{cov2,cov2}(1.,0.,0.,1.))

# g and g_inv must be inverse of each other, a `@check` that is only active when
# Gridap runs in debug mode
@test_throws AssertionError parse_metric_data(g_diag, SymTensorValue{con2}(1.,0.,1.))
@test_throws AssertionError parse_metric_data(sdg_diag, SymTensorValue{con2}(1.,0.,1.), g_diag)
@test_throws AssertionError compute_metric_data(g_diag, SymTensorValue{con2}(1.,0.,1.))

# `LinearAlgebra.I` parses into a `MetricData` storing no tensor, on which every
# operation selects its Euclidean version
md_I = parse_metric_data(I)
@test isidentity(md_I)
@test !isidentity(md_lazy) && !isidentity(md_eager)
@test md_I.g === I && md_I.g_inv === I && md_I.signed_sqrt_g === true
@test compute_metric_data(I) === md_I
@test compute_metric_data(md_I) === md_I
@test hodge_star(dx1, I)          === flat_hodge_star(dx1)
@test hodge_star(dx1, md_I)       === flat_hodge_star(dx1)
@test flat(v_gen, I)              === flat_flat(v_gen)
@test flat(e1, I)                 === flat_flat(e1)
@test flat(e1, md_I)              === flat_flat(e1)
@test sharp(omega_gen, I)         === flat_sharp(omega_gen)
@test sharp(omega_gen, md_I)      === flat_sharp(omega_gen)
@test alt_inner(dx1, dx1, I)      === flat_alt_inner(dx1, dx1)
@test alt_inner(dx1, dx1, md_I)   === flat_alt_inner(dx1, dx1)
@test hodge_star(TensorValue{cov2,cov2}(0.,-1.,1.,0.), I) === flat_hodge_star(TensorValue{cov2,cov2}(0.,-1.,1.,0.))


# ── Indefinite metric ─────────────────────────────────────────────────────────
#
# g = diag(-1,1) has det g = -1, so the volume density is 1 and the index that
# `signed_sqrt_g` carries makes it -1.

g_lor     = SymTensorValue{cov2,Float64,3}(-1.0, 0.0, 1.0)
g_inv_lor = inv(g_lor) # con2

@test compute_metric_data(g_lor).signed_sqrt_g         ≈ -1.0
@test compute_metric_data(g_inv_lor).signed_sqrt_g     ≈ -1.0
@test parse_metric_data(g_lor).signed_sqrt_g           ≈ -1.0
@test parse_metric_data(g_inv_lor).signed_sqrt_g       ≈ -1.0
@test compute_metric_data(g_lor, -1.0).signed_sqrt_g === -1.0

# A density of the wrong sign contradicts the index, a `@check` that is only
# active when Gridap runs in debug mode
@test_throws AssertionError parse_metric_data(g_lor, 1.0)
@test_throws AssertionError compute_metric_data(g_inv_lor, 1.0)

# `hodge_star` reads the density `abs(signed_sqrt_g)`, here 1: ⋆ of a 0-form is
# multiplication by it, and ⋆ of the top form by it times det g⁻¹
md_lor = compute_metric_data(g_lor)
@test collect(hodge_star(f0, md_lor).data)        ≈ [ 5.0]  atol=1e-12
@test collect(hodge_star(vol_coord, md_lor).data) ≈ [-1.0]  atol=1e-12

# ⋆ squares to sgn(det g)(-1)^{K(D-K)}, which is +1 on the 1-forms and -1 on
# the 0- and 2-forms of this metric
star2(ω) = hodge_star(hodge_star(ω, md_lor), md_lor)
@test collect(star2(f0).data)        ≈ -collect(f0.data)         atol=1e-12
@test collect(star2(vol_coord).data) ≈ -collect(vol_coord.data)  atol=1e-12
for ω in (dx1, dx2)
  @test collect(star2(ω).data) ≈ collect(ω.data)  atol=1e-12
end


# ── Inverse Hodge star ────────────────────────────────────────────────────────
#
# `⋆⁻¹ = sgn(det g)(-1)^{K(D-K)} ⋆` inverts `⋆` on either side, on a form and on
# a K-vector, whatever the index of the metric.

g2_gen = SymTensorValue{cov2}(2.0, 0.5, 3.0)
g3_gen = SymTensorValue{cov3}(2.0, 0.3, 0.1, 3.0, 0.2, 4.0)
g3_lor = SymTensorValue{cov3}(-1.0, 0.0, 0.0, 1.0, 0.0, 1.0)

for (D, mds) in ((2, (compute_metric_data(g2_gen), compute_metric_data(g_lor))),
                 (3, (compute_metric_data(g3_gen), compute_metric_data(g3_lor))))
  for K in 0:D, (Alias, Var) in ((ExteriorFormValue, Cov), (ExteriorVectorValue, Con))
    ω = Alias{K,Var(D)}(ntuple(i -> 1.0 + i, binomial(D,K)))

    for md in mds
      σ = sign(det(md.g)) * (isodd(K*(D-K)) ? -1 : 1)
      @test collect(inverse_hodge_star(hodge_star(ω, md), md).data) ≈ collect(ω.data)  atol=1e-12
      @test collect(hodge_star(inverse_hodge_star(ω, md), md).data) ≈ collect(ω.data)  atol=1e-12
      @test collect(inverse_hodge_star(ω, md).data) ≈ σ * collect(hodge_star(ω, md).data)  atol=1e-12
    end

    # The identity metric selects the Euclidean version, whose index is +1
    @test inverse_hodge_star(ω, I) === flat_inverse_hodge_star(ω)
    @test flat_inverse_hodge_star(flat_hodge_star(ω)) === ω
    @test flat_hodge_star(flat_inverse_hodge_star(ω)) === ω
  end
end

# A not-alternating tensor is inverse starred through its alternating part, and
# missing metric information is an error, not a silent Euclidean fallback
@test collect(inverse_hodge_star(TensorValue{cov2,cov2}(0.,-1.,1.,0.), md_eager).data) ≈
      collect(inverse_hodge_star(ExteriorFormValue{2,cov2}((1.,)), md_eager).data)  atol=1e-12
@test_throws ArgumentError inverse_hodge_star(dx1)
@test_throws DimensionMismatch inverse_hodge_star(dx1, md3)

# The aliases are the functions themselves
@test (⋆⁻¹) === inverse_hodge_star && (⋆ᴱ⁻¹) === flat_inverse_hodge_star
@test ⋆⁻¹(dx1, md_eager) === inverse_hodge_star(dx1, md_eager)
@test ⋆ᴱ⁻¹(dx1)         === flat_inverse_hodge_star(dx1)

end # module
