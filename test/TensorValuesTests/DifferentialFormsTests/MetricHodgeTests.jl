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
  SymTensorValue{2,Float64,3}(c*sa, -c*a*b, c*sb)
end

function _csphere_inv_metric(r, x)
  a, b   = tan(x[1]), tan(x[2])
  sa, sb = 1 + a^2, 1 + b^2
  ρ2     = 1 + a^2 + b^2
  c      = ρ2 / (r^2 * sa * sb)
  SymTensorValue{2,Float64,3}(c*sb, c*a*b, c*sa)
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
g_diag     = SymTensorValue{2,Float64,3}(p, 0.0, q)
g_inv_diag = SymTensorValue{2,Float64,3}(1/p, 0.0, 1/q)
sdg_diag   = sqrt(p*q)

dx1 = ExteriorFormValue{1,2}((1.0, 0.0))
dx2 = ExteriorFormValue{1,2}((0.0, 1.0))

star_dx1 = hodge_star(dx1, g_inv_diag, sdg_diag)
star_dx2 = hodge_star(dx2, g_inv_diag, sdg_diag)

@test typeof(star_dx1) <: ExteriorFormValue{1,2}
@test collect(star_dx1.data) ≈ [0.0, sqrt(q/p)]    atol=1e-12   # √(q/p) dx²
@test collect(star_dx2.data) ≈ [-sqrt(p/q), 0.0]   atol=1e-12   # -√(p/q) dx¹

# ⋆_g ∘ ⋆_g = (-1)^{K(D-K)} = -1 for K=1, D=2
@test collect(hodge_star(star_dx1, g_inv_diag, sdg_diag).data) ≈ collect((-1.0 * dx1).data)  atol=1e-12
@test collect(hodge_star(star_dx2, g_inv_diag, sdg_diag).data) ≈ collect((-1.0 * dx2).data)  atol=1e-12

# 0-form: ⋆_g f = f √det(g) vol_coord
f0    = ExteriorFormValue{0,2}((5.0,))
sf0_d = hodge_star(f0, g_inv_diag, sdg_diag)
@test typeof(sf0_d) <: ExteriorFormValue{2,2}
@test collect(sf0_d.data) ≈ [5.0 * sdg_diag]   atol=1e-12

# 2-form: ⋆_g vol_coord = 1/√det(g)
vol_coord = ExteriorFormValue{2,2}((1.0,))
svol_d    = hodge_star(vol_coord, g_inv_diag, sdg_diag)
@test typeof(svol_d) <: ExteriorFormValue{0,2}
@test collect(svol_d.data) ≈ [1.0/sdg_diag]   atol=1e-12

# ── Part 2: CubedSphere metric at a generic point ─────────────────────────────

# ⋆_g ∘ ⋆_g = -1 for K=1, D=2
for ω in [dx1, dx2]
  sω  = hodge_star(ω,  g_inv, sdg)
  ssω = hodge_star(sω, g_inv, sdg)
  @test collect(ssω.data) ≈ collect((-1.0 * ω).data)   atol=1e-10
end

# Generic 1-form
omega_gen = ExteriorFormValue{1,2}((0.7, -0.3))
s1 = hodge_star(omega_gen, g_inv, sdg)
s2 = hodge_star(s1,        g_inv, sdg)
@test collect(s2.data) ≈ collect((-1.0 * omega_gen).data)   atol=1e-10

# 0-form ⋆_g: result is a 2-form with coefficient f·√det(g)
f0_cs  = ExteriorFormValue{0,2}((3.0,))
sf0_cs = hodge_star(f0_cs, g_inv, sdg)
@test typeof(sf0_cs) <: ExteriorFormValue{2,2}
@test collect(sf0_cs.data) ≈ [3.0 * sdg]   atol=1e-10

# ── Part 3: Musical operators ─────────────────────────────────────────────────

e1 = VectorValue(1.0, 0.0)
e2 = VectorValue(0.0, 1.0)

# Flat metric: flat/sharp are identities
g_flat     = SymTensorValue{2,Float64,3}(1.0, 0.0, 1.0)
g_inv_flat = SymTensorValue{2,Float64,3}(1.0, 0.0, 1.0)

@test flat(e1, g_flat).data == (1.0, 0.0)
@test flat(e2, g_flat).data == (0.0, 1.0)
@test sharp(dx1, g_inv_flat) == e1
@test sharp(dx2, g_inv_flat) == e2

# Diagonal metric: flat(eᵢ, g) lowers with g
@test collect(flat(e1, g_diag).data) ≈ [p, 0.0]
@test collect(flat(e2, g_diag).data) ≈ [0.0, q]

# Round-trip: sharp(flat(v, g), g_inv) = v
for v in [e1, e2, VectorValue(1.5, -0.7)]
  v_rt = sharp(flat(v, g_diag), g_inv_diag)
  @test collect(v_rt.data) ≈ collect(v.data)   atol=1e-12
end

# Round-trip: flat(sharp(ω, g_inv), g) = ω
for ω in [dx1, dx2, ExteriorFormValue{1,2}((0.4, -0.9))]
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

end # module
