module MetricOperationsTests

using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Fields
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.TensorValues
using FillArrays
using Test
import LinearAlgebra

# The metric dependent operators on cell fields: spelled without a metric they
# operate with the metric of the triangulation they are evaluated on, and a
# metric given explicitly overrides it.

# ── A chart model of constant metric ──────────────────────────────────────────

# g = [[2, 1/2], [1/2, 1]], covariant as a metric must be. It is constant, so
# the operators have a closed form at every point.
const g = SymTensorValue{Cov(2)}(2.0,0.5,1.0)

struct ConstantMetric <: Field end
Gridap.Arrays.evaluate!(cache, ::ConstantMetric, x::Point{2}) = g

struct ChartModel{Dc,Dp,M,A} <: DiscreteModel{Dc,Dp}
  model::M
  cell_metric::A
  function ChartModel(model::DiscreteModel{Dc,Dp},cell_metric) where {Dc,Dp}
    new{Dc,Dp,typeof(model),typeof(cell_metric)}(model,cell_metric)
  end
end

Gridap.Geometry.get_grid(m::ChartModel)          = get_grid(m.model)
Gridap.Geometry.get_grid_topology(m::ChartModel) = get_grid_topology(m.model)
Gridap.Geometry.get_face_labeling(m::ChartModel) = get_face_labeling(m.model)
Gridap.Geometry.get_cell_metric(m::ChartModel)   = m.cell_metric

# The values a cell field takes at the first cell point of each of its cells.
function first_point_values(cf)
  a = cf(get_cell_points(get_triangulation(cf)))
  c = array_cache(a)
  [first(getindex!(c,a,cell)) for cell in 1:length(a)]
end

# The operator the operation tree holds, and the fields it is applied to.
operation_field(cf) = testitem(get_data(cf))

model  = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(2,2))
cmodel = ChartModel(model,Fill(ConstantMetric(),num_cells(model)))
Ω  = Triangulation(model)
Ωc = Triangulation(cmodel)

ωfun(x) = ExteriorFormValue{1,Cov(2)}(1.0 + 3x[1]*x[2], 2.0 + x[1]^2)
ηfun(x) = ExteriorFormValue{1,Cov(2)}(0.5, 1.0)
vfun(x) = VectorValue(1.0 + x[1], 2.0 + x[2])

ω, η, v    = map(f -> CellField(f,Ω ,PhysicalDomain()), (ωfun,ηfun,vfun))
ωc, ηc, vc = map(f -> CellField(f,Ωc,PhysicalDomain()), (ωfun,ηfun,vfun))

# ── A stock model operates with the Euclidean operator ────────────────────────

@test first_point_values(hodge_star(ω)) == first_point_values(flat_hodge_star(ω))
@test first_point_values(sharp(ω))      == first_point_values(flat_sharp(ω))
@test first_point_values(flat(v))       == first_point_values(flat_flat(v))
@test first_point_values(inner(ω,η))    == first_point_values(flat_inner(ω,η))
@test first_point_values(alt_inner(ω,η)) == first_point_values(flat_alt_inner(ω,η))

# The identity metric is not an operand: the Euclidean operator is applied to the
# fields alone, so a model declaring no metric evaluates no metric.
@test operation_field(hodge_star(ω)) isa Fields.OperationField{typeof(flat_hodge_star)}
@test operation_field(sharp(ω))      isa Fields.OperationField{typeof(flat_sharp)}
@test operation_field(flat(v))       isa Fields.OperationField{typeof(flat_flat)}
@test operation_field(inner(ω,η))    isa Fields.OperationField{typeof(flat_inner)}
@test operation_field(alt_inner(ω,η)) isa Fields.OperationField{typeof(flat_alt_inner)}
@test length(operation_field(inner(ω,η)).fields) == 2

# ── A chart model operates with the metric its model carries ──────────────────

md = compute_metric_data(g)

# The metric is an operand of the pointwise operator, after the fields it is
# spelled with, and only the data the operator reads is passed: `g` and `g_inv`
# for an inner product, `g_inv` alone to raise an index, and the volume factor
# only where it is part of the operation.
@test operation_field(inner(ωc,ηc)) isa Fields.OperationField{typeof(inner)}
@test length(operation_field(inner(ωc,ηc)).fields)     == 4
@test length(operation_field(alt_inner(ωc,ηc)).fields) == 4
@test length(operation_field(hodge_star(ωc)).fields)   == 4
@test length(operation_field(sharp(ωc)).fields)        == 2
@test length(operation_field(flat(vc)).fields)         == 2

@test all(first_point_values(hodge_star(ωc)) .≈ map(w -> hodge_star(w,md), first_point_values(ωc)))
@test all(first_point_values(sharp(ωc))      .≈ map(w -> sharp(w,md)     , first_point_values(ωc)))
@test all(first_point_values(flat(vc))       .≈ map(u -> flat(u,md)      , first_point_values(vc)))
@test all(first_point_values(inner(ωc,ηc))   .≈
          map((w,e) -> inner(w,e,md)    , first_point_values(ωc), first_point_values(ηc)))
@test all(first_point_values(alt_inner(ωc,ηc)) .≈
          map((w,e) -> alt_inner(w,e,md), first_point_values(ωc), first_point_values(ηc)))

# The metric reaches the quadrature points: with this g and η, ⟨ω,η⟩_g is the
# second component of ω, whose integral over the unit square is 2 + 1/3.
# `Measure` integrates against the euclidean volume element, so the integral
# carries no √|det g| factor.
@test sum(∫( inner(ωc,ηc) )Measure(Ωc,2)) ≈ 7/3

# A face of a chart model needs the metric induced on it, which the model does
# not carry, so the operators on it fail rather than read the cell metric.
Γc = BoundaryTriangulation(cmodel)
ωγ = CellField(ωfun,Γc,PhysicalDomain())
ηγ = CellField(ηfun,Γc,PhysicalDomain())
@test_throws ErrorException hodge_star(ωγ)
@test_throws ErrorException inner(ωγ,ηγ)

# ── An explicitly given metric overrides the triangulation's ──────────────────

g_cf     = CellField(x -> g     , Ωc, PhysicalDomain())
g_inv_cf = CellField(x -> inv(g), Ωc, PhysicalDomain())

# `LinearAlgebra.I` asks for the Euclidean operator on any model
@test first_point_values(hodge_star(ωc,LinearAlgebra.I)) == first_point_values(flat_hodge_star(ωc))
@test first_point_values(sharp(ωc,LinearAlgebra.I))      == first_point_values(flat_sharp(ωc))
@test first_point_values(flat(vc,LinearAlgebra.I))       == first_point_values(flat_flat(vc))
@test first_point_values(inner(ωc,ηc,LinearAlgebra.I))   == first_point_values(flat_inner(ωc,ηc))
@test first_point_values(alt_inner(ωc,ηc,LinearAlgebra.I)) == first_point_values(flat_alt_inner(ωc,ηc))

# a cell field of metric tensors, a constant tensor and the inverse pair all
# reproduce the implicit metric of this chart
@test first_point_values(hodge_star(ωc,g_cf))            ≈ first_point_values(hodge_star(ωc))
@test first_point_values(hodge_star(ωc,g))               ≈ first_point_values(hodge_star(ωc))
@test first_point_values(hodge_star(ωc,g_cf,g_inv_cf))   ≈ first_point_values(hodge_star(ωc))
@test first_point_values(inner(ωc,ηc,g_cf))              ≈ first_point_values(inner(ωc,ηc))
@test first_point_values(alt_inner(ωc,ηc,g))             ≈ first_point_values(alt_inner(ωc,ηc))
@test first_point_values(sharp(ωc,g_inv_cf))             ≈ first_point_values(sharp(ωc))

# a metric constant per cell is given as the array of its cell values
@test first_point_values(sharp(ωc,Fill(inv(g),num_cells(cmodel)))) ≈ first_point_values(sharp(ωc))

# ── The tensor inner product of fields without variance is unchanged ──────────

us = CellField(x -> x[1]*x[2], Ω, PhysicalDomain())
ut = CellField(x -> TensorValue(1.0 + x[1],0.0,0.0,x[2]), Ω, PhysicalDomain())
vt = CellField(x -> TensorValue(1.0,0.0,0.0,2.0)        , Ω, PhysicalDomain())
tfun(x) = TensorValue(1.0,0.0,0.0,2.0)

@test first_point_values(inner(ut,vt)) == first_point_values(flat_inner(ut,vt))
@test first_point_values(inner(ut,tfun)) == first_point_values(inner(ut,vt))
@test first_point_values(inner(tfun,ut)) == first_point_values(inner(ut,vt))
@test first_point_values(inner(us,2.0)) == 2 .* first_point_values(us)
@test first_point_values(inner(2.0,us)) == 2 .* first_point_values(us)

# The Leibniz rule holds of the euclidean inner product as of the varianced one:
# ⟨ut,vt⟩ = 1 + x¹ + 2x², whose partial derivatives are (1,2).
@test all(d -> (d[1],d[2]) == (1.0,2.0), first_point_values(∂(inner(ut,vt))))

# It holds of the lazy arrays an FE function and an FE basis carry, which take
# their own rules, as well as of a single field.
W  = FESpace(model,ReferenceFE(lagrangian,TensorValue{2,2,Float64,4},1))
uh = interpolate(x -> TensorValue(1.0 + x[1],0.0,0.0,x[2]), W)
@test all(d -> d[1] ≈ 1.0 && d[2] ≈ 2.0, first_point_values(∂(inner(uh,vt))))
@test all(d -> d isa MultiValue{Tuple{Cov(2)}}, first_point_values(∂(inner(get_fe_basis(W),vt))))

# Differentiating through an injected metric would need the derivatives of the
# metric, which the Leibniz rule of a product of two fields does not take.
utc = CellField(x -> SymTensorValue{Cov(2)}(1.0 + x[1],0.0,x[2]), Ωc, PhysicalDomain())
@test_throws ErrorException ∂(inner(utc,utc))

# ── FE assembly on a stock model is untouched ─────────────────────────────────

V  = FESpace(model,ReferenceFE(lagrangian,VectorValue{2,Float64},1),dirichlet_tags="boundary")
dΩ = Measure(Ω,2)
a(u,w) = ∫( inner(∇(u),∇(w)) )dΩ
l(w)   = ∫( inner(VectorValue(1.0,1.0),w) )dΩ
uh = solve(AffineFEOperator(a,l,V,V))
@test sum(∫( inner(uh,uh) )dΩ) > 0

end # module
