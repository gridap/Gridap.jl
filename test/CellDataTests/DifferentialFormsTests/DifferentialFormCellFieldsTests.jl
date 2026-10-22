module DifferentialFormCellFieldsTests

using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: GenericField
using Test
import LinearAlgebra
using Gridap.CellData: GenericCellField, get_data
using FillArrays: Fill

# Form-valued CellFields: construction, evaluation, form-level and
# value-level operators, the flat_vol_coeff/to_1form/from_1form bridges, and
# near-zero-allocation cell traversal.
#
# The 1-form ω = 3xy dx¹ + x² dx²  is used throughout.
# Its exterior derivative: dω = (∂x²/∂x − ∂3xy/∂y) dx¹∧dx² = −x dx¹∧dx²
# Codifferential (sign = -1 for D=2, K=1):
#   ⋆(ω₁dx¹+ω₂dx²) = (−ω₂)dx¹ + ω₁dx²  (convention: ⋆dx¹=dx², ⋆dx²=−dx¹)
#   ⋆ω = (−x²)dx¹ + (3xy)dx²
#   d⋆ω = (∂(3xy)/∂x − ∂(−x²)/∂y) dx¹∧dx² = 3y dx¹∧dx²
#   δω = −1 × ⋆(3y vol) = −3y;  at (0.125,0.125): −0.375

# The values a cell field takes at the cell points of its own triangulation.
function cell_values(cf)
  a = cf(get_cell_points(get_triangulation(cf)))
  c = array_cache(a)
  [getindex!(c, a, i) for i in 1:length(a)]
end

# Element type of the values a form valued cell field produces, read off the
# values themselves rather than off a declared type parameter.
form_eltype(cf) = eltype(first(cell_values(cf)))

# ── Setup ─────────────────────────────────────────────────────────────────────

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0), (4,4))
Ω     = Triangulation(model)

# The form is one field returning an ExteriorFormValue, so the differential
# operators differentiate it as a whole and need no access to its components.
ωfun(x) = CoVectorValue(3*x[1]*x[2], x[1]^2)
ω = CellField(ωfun, Ω, PhysicalDomain())

@test get_triangulation(ω) === Ω
@test DomainStyle(ω) == PhysicalDomain()
@test form_eltype(ω) <: ExteriorFormValue{1,cov2}

# ── Evaluation at cell points ─────────────────────────────────────────────────

pts   = get_cell_points(Ω)
ωx    = ω(pts)                       # LazyArray — nothing computed yet

cω    = array_cache(ωx)
val0  = getindex!(cω, ωx, 1)         # first cell: ExteriorFormValue per quad pt
@test typeof(val0) <: AbstractArray{<:CoVectorValue{cov2}}

# ── exterior_derivative ───────────────────────────────────────────────────────

dω = exterior_derivative(ω)
@test form_eltype(dω) <: ExteriorFormValue{2,cov2}

dωx   = dω(pts)
cdω   = array_cache(dωx)
dval0 = getindex!(cdω, dωx, 1)
@test typeof(dval0) <: AbstractArray{<:ExteriorFormValue{2,cov2}}

# Verify value: dω = -x dx¹∧dx² → at (0.125,0.125): coefficient = -0.125
x_c = Point(0.125, 0.125)
ωfld = GenericField(ωfun)
dv = evaluate(exterior_derivative(ωfld), x_c)
@test dv.data[1] ≈ -0.125   atol=1e-12

# ── codifferential ────────────────────────────────────────────────────────────

δω = codifferential(ω)
@test form_eltype(δω) <: ExteriorFormValue{0,cov2}

δωx   = δω(pts)
cδω   = array_cache(δωx)
δval0 = getindex!(cδω, δωx, 1)
@test typeof(δval0) <: AbstractArray{<:ExteriorFormValue{0,cov2}}

# Verify value: δω = -3y → at (0.125,0.125): -0.375
δv = evaluate(codifferential(ωfld, LinearAlgebra.I), x_c)
@test δv.data[1] ≈ -0.375   atol=1e-8

# `δ` of a 0-form is the zero (-1)-form, which has no value type
ω0 = CellField(x -> ExteriorFormValue{0,cov2}(x[1]*x[2]), Ω, PhysicalDomain())
@test_throws ErrorException codifferential(ω0)(pts)

# ── hodge_star (value level via Operation) ────────────────────────────────────

star_ω = flat_hodge_star(ω)
@test form_eltype(star_ω) <: ExteriorFormValue{1,cov2}

star_valx   = star_ω(pts)
csv         = array_cache(star_valx)
sv0         = getindex!(csv, star_valx, 1)
@test typeof(sv0) <: AbstractArray{<:CoVectorValue{cov2}}

# Chaining: the star is a pointwise operation and the exterior derivative
# differentiates the starred field as a black box, so d(⋆ω) is a 2-form.
d_star_ω = exterior_derivative(star_ω)
@test form_eltype(d_star_ω) <: ExteriorFormValue{2,cov2}

# δω = (-1)^{D(K-1)+1} ⋆d⋆ω, with D=2 and K=1 the sign is -1
@test (-flat_hodge_star(evaluate(exterior_derivative(Operation(flat_hodge_star)(ωfld)), x_c))).data[1] ≈
      δv.data[1]   atol=1e-8

# ── wedge product ─────────────────────────────────────────────────────────────

# η = dx² (constant 1-form)
ηfun(x) = ExteriorFormValue{1,Cov(2)}((0.0, 1.0))
η = CellField(ηfun, Ω, PhysicalDomain())

wedge_cf = ω ∧ η
# ω ∧ η = (3xy dx¹ + x² dx²) ∧ dx² = 3xy dx¹∧dx²  (since dx²∧dx²=0)
wedge_x  = wedge_cf(pts)
cw       = array_cache(wedge_x)
wval0    = getindex!(cw, wedge_x, 1)
@test typeof(wval0) <: AbstractArray{<:ExteriorFormValue{2,cov2}}

# ── Arithmetic ────────────────────────────────────────────────────────────────

two_omega = 2.0 * ω
two_omega_x = two_omega(pts)
ct          = array_cache(two_omega_x)
tval0       = getindex!(ct, two_omega_x, 1)
@test tval0 ≈ 2.0 .* val0

# ── flat and sharp (value-level via Operation) ────────────────────────────────

g_cf     = CellField(x -> SymTensorValue{cov2,Float64,3}(1.0,0.0,1.0), Ω, PhysicalDomain())
g_inv_cf = CellField(x -> SymTensorValue{con2,Float64,3}(1.0,0.0,1.0), Ω, PhysicalDomain())
v_cf     = CellField(x -> VectorValue(x[1], x[2]), Ω, PhysicalDomain())

flat_cf  = flat(v_cf, g_cf)
flatx    = flat_cf(pts)
cflat    = array_cache(flatx)
fval0    = getindex!(cflat, flatx, 1)
@test typeof(fval0) <: AbstractArray   # CoVectorValue{2} elements

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
@test typeof(ipv0) <: AbstractArray   # ExteriorFormValue{0,2} elements

# ── koszul ────────────────────────────────────────────────────────────────────
# κ_x(ω) = ι_x(ω) needs the evaluation point, which the identity field supplies

ω_id_fun(x) = CoVectorValue(x[1], x[2])
ω_id_cf = CellField(ω_id_fun, Ω, PhysicalDomain())
κω_cf = koszul(ω_id_cf)
@test form_eltype(κω_cf) <: ExteriorFormValue{0,cov2}

κωx   = κω_cf(pts)
cκ    = array_cache(κωx)
κval0 = getindex!(cκ, κωx, 1)
@test typeof(κval0) <: AbstractArray{<:ExteriorFormValue{0,cov2}}

# κ(x¹dx¹ + x²dx²) = (x¹)² + (x²)², the squared norm of the position. The cell
# field is in the physical domain, so that is where the identity field reads x.
qp = get_data(change_domain(pts, PhysicalDomain()))[1]
@test all(κval0[i].data[1] ≈ qp[i]⋅qp[i] for i in eachindex(κval0, qp))

# ── flat_vol_coeff (CellField level) ──────────────────────────────────────────

ω2_cf = CellField(x -> ExteriorFormValue{2,Cov(2)}((7.0,)), Ω, PhysicalDomain())
scalar_cf = flat_vol_coeff(ω2_cf)
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
@test typeof(ωvval0) <: AbstractArray{<:CoVectorValue{cov2}}

v_back = from_1form(ω_id_cf)
@test v_back isa CellField

vx    = v_back(pts)
cv    = array_cache(vx)
vval0 = getindex!(cv, vx, 1)
@test typeof(vval0) <: AbstractArray{<:VectorValue{2}}

# ── Integration: ∫(flat_vol_coeff(ω ∧ ⋆ω)) dΩ equals ∫‖ω‖² dΩ ───────────────
# ω = x¹ dx¹ + x² dx² on [0,1]²
# ‖ω‖² = x¹² + x²²; ∫₀¹∫₀¹(x²+y²) dx dy = 1/3+1/3 = 2/3
dΩ = Measure(Ω, 4)
norm2 = sum(∫(flat_vol_coeff(ω_id_cf ∧ flat_hodge_star(ω_id_cf))) * dΩ)
@test norm2 ≈ 2/3   atol=1e-4

# ── Near-zero-allocation cell traversal ──────────────────────────────────────
# Build a ReferenceDomain form cell field directly to avoid the
# change_domain(Physical→Reference) step; this checks that the exterior
# derivative evaluation itself is (nearly) allocation-free.

dω_ref = GenericCellField(
  lazy_map(Broadcasting(exterior_derivative), get_data(ω)),
  Ω, ReferenceDomain())

dωx_ref = dω_ref(pts)                  # same domain — no change_domain step
cdω_ref = array_cache(dωx_ref)
getindex!(cdω_ref, dωx_ref, 1)          # warm up

allocs = @allocated for i in 1:num_cells(Ω)
  getindex!(cdω_ref, dωx_ref, i)
end
@test allocs / num_cells(Ω) < 500

# ── Derivatives of linear combinations: FE functions and non affine maps ─────

V  = TestFESpace(model, ReferenceFE(lagrangian,Float64,2))
uh = FEFunction(V, rand(num_free_dofs(V)))
du = get_fe_basis(V)

for fe_cf in (uh, du)
  dfe  = exterior_derivative(fe_cf)(pts)
  cdfe = array_cache(dfe)
  @test eltype(getindex!(cdfe, dfe, 1)) <: ExteriorFormValue{1,cov2}

  # The exterior derivative of a 0-form is its gradient read as a 1-form.
  @test all(all(a .≈ b) for (a,b) in
            zip(cell_values(𝑑(fe_cf)), cell_values(to_1form(∇(fe_cf)))))
end

# The jacobian of a bilinear cell map holds the same numbers as `∇` and has
# opposite variances, so `pinvJt` accepts it.
Ωu  = Triangulation(UnstructuredDiscreteModel(model))
cmu = get_cell_map(Ωu)
q   = Point(0.3,0.7)
Jt  = evaluate(lazy_map(partial_derivatives, cmu)[1], q)
@test isa(Jt, MultiValue{Tuple{cov2,con2}})
@test get_array(Jt) ≈ get_array(evaluate(∇(cmu[1]), q))
@test get_array(pinvJt(Jt)) ≈ get_array(pinvJt(evaluate(∇(cmu[1]), q)))

# ── Derivatives in the reference domain ──────────────────────────────────────
# A reference domain cell field takes reference coordinates and returns values
# in the physical basis, so differentiating it transports its derivation
# indices alone, which `pushforward_partials` does, as `∇` does with `push_∇`.
# The carrier of the result is the one `derivatives_type` prescribes, in either
# domain.

∂ω = partial_derivatives(ω)
@test DomainStyle(∂ω) == PhysicalDomain()
@test form_eltype(∂ω) === derivatives_type(CoVectorValue{cov2,Float64}, x_c)

@test DomainStyle(uh) == ReferenceDomain()
∂uh = partial_derivatives(uh)
@test DomainStyle(∂uh) == ReferenceDomain()
@test form_eltype(∂uh) === derivatives_type(Float64, x_c)
@test all(all(get_array(ai) ≈ get_array(bi) for (ai,bi) in zip(as,bs))
          for (as,bs) in zip(cell_values(∂uh), cell_values(∇(uh))))

# Pushing a higher order derivative forward needs the derivatives of the cell
# map itself, as `push_∇∇` does.
@test_throws ErrorException second_partial_derivatives(uh)

# The operators agree with their physical twins, on a mesh whose cell maps are
# not affine as well as on the affine one.
model_na = CartesianDiscreteModel((0.0,1.0,0.0,1.0), (4,4),
                                  map = x -> Point(x[1]+0.3*x[1]*x[2], x[2]+0.1*x[1]))

# The Lie derivative reads the partial derivatives of both its arguments and
# needs no metric, so the generic `Fields` operation already carries it.
uflow(x) = ContraVectorValue(x[1] + 2*x[2], x[2]^2)
lie_along(a) = 𝓛(CellField(uflow, get_triangulation(a), PhysicalDomain()), a)

for Ωi in (Ω, Triangulation(model_na))
  ω_phys = CellField(ωfun, Ωi, PhysicalDomain())
  ω_ref  = change_domain(ω_phys, Ωi, ReferenceDomain())
  @test DomainStyle(ω_ref) == ReferenceDomain()

  for op in (exterior_derivative, codifferential, koszul, partial_derivatives, lie_along)
    @test form_eltype(op(ω_ref)) === form_eltype(op(ω_phys))
    @test all(all(a .≈ b) for (a,b) in
              zip(cell_values(op(ω_phys)), cell_values(op(ω_ref))))
  end
end

# ── The FEEC spelling of a bilinear form ─────────────────────────────────────
# On a 0-form space `∫ 𝑑u ∧ ⋆𝑑v` is the stiffness matrix `∫ ∇u⋅∇v`, which pins
# the pushforward: the FE basis lives in the reference domain.

U = TrialFESpace(V)
A_dg = assemble_matrix((u,v) -> ∫( flat_vol_coeff(𝑑(u) ∧ flat_hodge_star(𝑑(v))) )dΩ, U, V)
A_eu = assemble_matrix((u,v) -> ∫( ∇(u)⋅∇(v) )dΩ, U, V)
@test A_dg ≈ A_eu

# On a 0-form Cartan's formula is `𝓛ᵤf = ι_u(𝑑f)`, which pins 𝓛 on a FE basis
# living in the reference domain.
uf_cf = CellField(uflow, Ω, PhysicalDomain())
B_lie = assemble_matrix((u,v) -> ∫( 𝓛(uf_cf,u)*v )dΩ, U, V)
B_cartan = assemble_matrix((u,v) -> ∫( interior_product(uf_cf, 𝑑(u))*v )dΩ, U, V)
@test B_lie ≈ B_cartan

end # module
