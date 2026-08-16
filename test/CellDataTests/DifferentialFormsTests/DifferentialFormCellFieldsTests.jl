module DifferentialFormCellFieldsTests
# DifferentialFormCellField: construction, evaluation, form-level and
# value-level operators, the vol_coeff/to_1form/from_1form bridges, and
# near-zero-allocation cell traversal.
#
# The 1-form ω = 3xy dx¹ + x² dx²  is used throughout.
# Its exterior derivative: dω = (∂x²/∂x − ∂3xy/∂y) dx¹∧dx² = −x dx¹∧dx²
# Codifferential (sign = -1 for D=2, K=1):
#   ⋆(ω₁dx¹+ω₂dx²) = (−ω₂)dx¹ + ω₁dx²  (convention: ⋆dx¹=dx², ⋆dx²=−dx¹)
#   ⋆ω = (−x²)dx¹ + (3xy)dx²
#   d⋆ω = (∂(3xy)/∂x − ∂(−x²)/∂y) dx¹∧dx² = 3y dx¹∧dx²
#   δω = −1 × ⋆(3y vol) = −3y;  at (0.125,0.125): −0.375

using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: GenericField
using Test

# ── Setup ─────────────────────────────────────────────────────────────────────

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0), (4,4))
Ω     = Triangulation(model)

# Build ω as a DifferentialForm so its components are differentiable
ω_form = DifferentialForm{1,2}((
  GenericField(x -> 3*x[1]*x[2]),   # ω₁ = 3xy
  GenericField(x -> x[1]^2)          # ω₂ = x²
))
ω = DifferentialFormCellField{1,2}(ω_form, Ω, PhysicalDomain())

@test get_triangulation(ω) === Ω
@test DomainStyle(ω) == PhysicalDomain()
@test typeof(ω) <: DifferentialFormCellField{1,2}

# ── Evaluation at cell points ─────────────────────────────────────────────────

pts   = get_cell_points(Ω)
ωx    = ω(pts)                       # LazyArray — nothing computed yet

cω    = array_cache(ωx)
val0  = getindex!(cω, ωx, 1)         # first cell: DifferentialFormValue per quad pt
@test typeof(val0) <: AbstractArray{<:DifferentialFormValue{1,2}}

# ── exterior_derivative ───────────────────────────────────────────────────────

dω = exterior_derivative(ω)
@test typeof(dω) <: DifferentialFormCellField{2,2}

dωx   = dω(pts)
cdω   = array_cache(dωx)
dval0 = getindex!(cdω, dωx, 1)
@test typeof(dval0) <: AbstractArray{<:DifferentialFormValue{2,2}}

# Verify value: dω = -x dx¹∧dx² → at (0.125,0.125): coefficient = -0.125
x_c = Point(0.125, 0.125)
dω_fld = ExteriorDerivativeForm(ω_form)
cache_d = return_cache(dω_fld, x_c)
dv = evaluate!(cache_d, dω_fld, x_c)
@test dv.data[1] ≈ -0.125   atol=1e-12

# ── codifferential ────────────────────────────────────────────────────────────

δω = codifferential(ω)
@test typeof(δω) <: DifferentialFormCellField{0,2}

δωx   = δω(pts)
cδω   = array_cache(δωx)
δval0 = getindex!(cδω, δωx, 1)
@test typeof(δval0) <: AbstractArray{<:DifferentialFormValue{0,2}}

# Verify value: δω = -3y → at (0.125,0.125): -0.375
δω_fld = codifferential(ω_form)
cache_δ = return_cache(δω_fld, x_c)
δv = evaluate!(cache_δ, δω_fld, x_c)
@test δv.data[1] ≈ -0.375   atol=1e-8

# ── hodge_star_form ───────────────────────────────────────────────────────────

star_ω = hodge_star_form(ω)
@test typeof(star_ω) <: DifferentialFormCellField{1,2}

# Chaining: d(⋆ω) should give back a 2-form CellField
d_star_ω = exterior_derivative(star_ω)
@test typeof(d_star_ω) <: DifferentialFormCellField{2,2}

# ── hodge_star (value level via Operation) ────────────────────────────────────

star_val_cf = hodge_star(ω)
star_valx   = star_val_cf(pts)
csv         = array_cache(star_valx)
sv0         = getindex!(csv, star_valx, 1)
@test typeof(sv0) <: AbstractArray   # DifferentialFormValue{1,2} elements

# ── wedge product ─────────────────────────────────────────────────────────────

# η = dx² (constant 1-form as a uniform DifferentialForm)
η_form = DifferentialForm{1,2}((
  GenericField(x -> 0.0),
  GenericField(x -> 1.0)
))
η = DifferentialFormCellField{1,2}(η_form, Ω, PhysicalDomain())

wedge_cf = ω ∧ η
# ω ∧ η = (3xy dx¹ + x² dx²) ∧ dx² = 3xy dx¹∧dx²  (since dx²∧dx²=0)
wedge_x  = wedge_cf(pts)
cw       = array_cache(wedge_x)
wval0    = getindex!(cw, wedge_x, 1)
@test typeof(wval0) <: AbstractArray   # DifferentialFormValue{2,2} elements

# ── Arithmetic ────────────────────────────────────────────────────────────────

two_omega = 2.0 * ω
two_omega_x = two_omega(pts)
ct          = array_cache(two_omega_x)
tval0       = getindex!(ct, two_omega_x, 1)
@test typeof(tval0) <: AbstractArray

# ── flat and sharp (value-level via Operation) ────────────────────────────────

g_cf     = CellField(x -> SymTensorValue{2,Float64,3}(1.0,0.0,1.0), Ω, PhysicalDomain())
g_inv_cf = CellField(x -> SymTensorValue{2,Float64,3}(1.0,0.0,1.0), Ω, PhysicalDomain())
v_cf     = CellField(x -> VectorValue(x[1], x[2]), Ω, PhysicalDomain())

flat_cf  = flat(v_cf, g_cf)
flatx    = flat_cf(pts)
cflat    = array_cache(flatx)
fval0    = getindex!(cflat, flatx, 1)
@test typeof(fval0) <: AbstractArray   # DifferentialFormValue{1,2} elements

sharp_cf = sharp(ω, g_inv_cf)
sharpx   = sharp_cf(pts)
csharp   = array_cache(sharpx)
shval0   = getindex!(csharp, sharpx, 1)
@test typeof(shval0) <: AbstractArray  # VectorValue{2} elements

# ── interior product ──────────────────────────────────────────────────────────

ip_cf = interior_product(v_cf, ω)
ipx   = ip_cf(pts)
cip   = array_cache(ipx)
ipv0  = getindex!(cip, ipx, 1)
@test typeof(ipv0) <: AbstractArray   # DifferentialFormValue{0,2} elements

# ── koszul ────────────────────────────────────────────────────────────────────

ω_id_form = DifferentialForm{1,2}((GenericField(x -> x[1]), GenericField(x -> x[2])))
ω_id_cf   = DifferentialFormCellField{1,2}(ω_id_form, Ω, PhysicalDomain())
κω_cf = koszul(ω_id_cf)
@test typeof(κω_cf) <: DifferentialFormCellField{0,2}

κωx   = κω_cf(pts)
cκ    = array_cache(κωx)
κval0 = getindex!(cκ, κωx, 1)
@test typeof(κval0) <: AbstractArray{<:DifferentialFormValue{0,2}}

# ── vol_coeff (CellField level) ───────────────────────────────────────────────

ω2_form = DifferentialForm{2,2}((GenericField(x -> 7.0),))
ω2_cf   = DifferentialFormCellField{2,2}(ω2_form, Ω, PhysicalDomain())
scalar_cf = vol_coeff(ω2_cf)
@test scalar_cf isa CellField

svals  = scalar_cf(pts)
cs     = array_cache(svals)
sv1    = getindex!(cs, svals, 1)
@test all(v ≈ 7.0 for v in sv1)

# ── to_1form / from_1form (CellField level) ───────────────────────────────────

ω_v_cf = to_1form(v_cf)
@test ω_v_cf isa CellField

ωvx    = ω_v_cf(pts)
cωv    = array_cache(ωvx)
ωvval0 = getindex!(cωv, ωvx, 1)
@test typeof(ωvval0) <: AbstractArray{<:DifferentialFormValue{1,2}}

v_back = from_1form(ω_id_cf)
@test v_back isa CellField

vx    = v_back(pts)
cv    = array_cache(vx)
vval0 = getindex!(cv, vx, 1)
@test typeof(vval0) <: AbstractArray{<:VectorValue{2}}

# ── Integration: ∫(vol_coeff(ω ∧ ⋆ω)) dΩ equals ∫‖ω‖² dΩ ────────────────────
# ω = x¹ dx¹ + x² dx² on [0,1]²
# ‖ω‖² = x¹² + x²²; ∫₀¹∫₀¹(x²+y²) dx dy = 1/3+1/3 = 2/3
dΩ = Measure(Ω, 4)
norm2 = sum(∫(vol_coeff(ω_id_cf ∧ hodge_star(ω_id_cf))) * dΩ)
@test norm2 ≈ 2/3   atol=1e-4

# ── Near-zero-allocation cell traversal ──────────────────────────────────────
# Build a ReferenceDomain DifferentialFormCellField directly to avoid the
# change_domain(Physical→Reference) step; this checks that the lazy
# ExteriorDerivativeForm evaluation itself is (nearly) allocation-free.

dω_ref = DifferentialFormCellField{2,2}(
  lazy_map(Broadcasting(exterior_derivative), Gridap.CellData.get_data(ω)),
  Ω, ReferenceDomain())

dωx_ref = dω_ref(pts)                  # same domain — no change_domain step
cdω_ref = array_cache(dωx_ref)
getindex!(cdω_ref, dωx_ref, 1)          # warm up

allocs = @allocated for i in 1:num_cells(Ω)
  getindex!(cdω_ref, dωx_ref, i)
end
@test allocs / num_cells(Ω) < 500

end # module
