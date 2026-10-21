module CellMetricsTests

using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.TensorValues
using Test
import LinearAlgebra

# The metric of a model, as the arguments the metric aware operators take, on a
# stock Cartesian model and on a mock chart model.

# ── A chart model ─────────────────────────────────────────────────────────────

# The metric of the chart, g = [[c + x¹, 3/10], [3/10, 3/2]], covariant as a
# metric must be. The constant c is one per cell, so that a cell metric read at
# the wrong cell shows.
struct CovChartMetric <: Field
  c::Float64
end
Gridap.Arrays.evaluate!(cache, m::CovChartMetric, x::Point{2}) =
  SymTensorValue{Cov(2)}(m.c + x[1], 0.3, 1.5)

# The same metric with its dimension left un-varianced, which carries no
# variance for the operators to read.
struct ChartMetric <: Field
  c::Float64
end
Gridap.Arrays.evaluate!(cache, m::ChartMetric, x::Point{2}) =
  SymTensorValue{2}(m.c + x[1], 0.3, 1.5)

# A contravariant metric tensor, which is an inverse metric, not a metric.
struct ConChartMetric <: Field
  c::Float64
end
Gridap.Arrays.evaluate!(cache, m::ConChartMetric, x::Point{2}) =
  SymTensorValue{Con(2)}(m.c + x[1], 0.3, 1.5)

# The closed form of g⁻¹, with det g = (3/2)(c + x¹) - 9/100.
struct InvChartMetric <: Field
  c::Float64
end
function Gridap.Arrays.evaluate!(cache, m::InvChartMetric, x::Point{2})
  g11 = m.c + x[1]
  d = 1.5*g11 - 0.09
  SymTensorValue{Con(2)}(1.5/d, -0.3/d, g11/d)
end

struct ChartModel{Dc,Dp,M,A,B} <: DiscreteModel{Dc,Dp}
  model::M
  cell_metric::A
  cell_inv_metric::B
  function ChartModel(model::DiscreteModel{Dc,Dp},cell_metric,cell_inv_metric=nothing) where {Dc,Dp}
    new{Dc,Dp,typeof(model),typeof(cell_metric),typeof(cell_inv_metric)}(
      model,cell_metric,cell_inv_metric)
  end
end

Gridap.Geometry.get_grid(m::ChartModel)          = get_grid(m.model)
Gridap.Geometry.get_grid_topology(m::ChartModel) = get_grid_topology(m.model)
Gridap.Geometry.get_face_labeling(m::ChartModel) = get_face_labeling(m.model)
Gridap.Geometry.get_cell_metric(m::ChartModel)   = m.cell_metric

# Only a model holding a closed form for g⁻¹ answers for it; the others leave it
# to be derived from the metric.
Gridap.Geometry.get_cell_inv_metric(m::ChartModel{Dc,Dp,M,A,<:AbstractArray}) where {Dc,Dp,M,A} =
  m.cell_inv_metric

cell_metric_fields(metric_field,model) =
  [metric_field(Float64(cell)) for cell in 1:num_cells(model)]
chart_model(model,metric_field) = ChartModel(model,cell_metric_fields(metric_field,model))

# The values a cell field takes at the first cell point of each of its cells.
function first_point_values(cf)
  a = cf(get_cell_points(get_triangulation(cf)))
  c = array_cache(a)
  [first(getindex!(c,a,cell)) for cell in 1:length(a)]
end

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(2,2))
Ω     = Triangulation(model)
ω     = ExteriorFormValue{1,2}(1.0,2.0)

# ── A stock model has the identity metric ─────────────────────────────────────

@test get_cell_metric(model) === LinearAlgebra.I
@test get_cell_inv_metric(model) === LinearAlgebra.I

@test get_metric_data(Ω) === (LinearAlgebra.I,)
@test get_metric_data(BoundaryTriangulation(model)) === (LinearAlgebra.I,)
@test get_metric_data(SkeletonTriangulation(model)) === (LinearAlgebra.I,)

# The identity metric reproduces the metric free operators
@test hodge_star(ω,get_metric_data(Ω)...) == flat_hodge_star(ω)
@test alt_inner(ω,ω,get_metric_data(Ω)...) == flat_alt_inner(ω,ω)

# ── A chart model carries its own metric ──────────────────────────────────────

# The cells of this mesh are ordered x¹ first, so the first point of cells 1 to 4
# has x¹ = 0, 1/2, 0, 1/2, and the metric there is cell + x¹.
expected_g_11 = [1.0, 2.5, 3.0, 4.5]

cmodel = chart_model(model,CovChartMetric)
Ωc = Triangulation(cmodel)

@test length(get_cell_metric(cmodel)) == num_cells(cmodel)
@test length(get_cell_inv_metric(cmodel)) == num_cells(cmodel)

g, g_inv, sqrt_abs_g = get_metric_data(Ωc)

@test DomainStyle(g) == PhysicalDomain()
@test DomainStyle(g_inv) == PhysicalDomain()
@test DomainStyle(sqrt_abs_g) == PhysicalDomain()

gv, giv, sv = map(first_point_values, (g,g_inv,sqrt_abs_g))

@test all(v -> v isa SymTensorValue{Cov(2),Float64}, gv)
@test all(v -> v isa SymTensorValue{Con(2),Float64}, giv)
@test [v[1,1] for v in gv] == expected_g_11
@test all(cell -> get_array(giv[cell]) ≈ inv(get_array(gv[cell])), 1:num_cells(Ωc))
@test all(cell -> sv[cell] ≈ sqrt(abs(det(get_array(gv[cell])))), 1:num_cells(Ωc))

# The three are what the pointwise operators ask for, in that order. The
# inverse and the determinant are checked against each other by MetricData.
md = compute_metric_data(gv[1],giv[1],sv[1])
@test md.g == gv[1]
@test md.g_inv == giv[1]
@test md.sqrt_abs_g == sv[1]
@test hodge_star(ω,gv[1],giv[1],sv[1]) ≈ hodge_star(ω,md)

# The cell fields reach the quadrature points. `Measure` scales the weights of a
# scalar integrand by `sqrt_abs_g`, which dividing by it cancels, leaving
# ∫ det g = ∫ (3/2)(cell + x¹) - 9/100 over the unit square
dΩc = Measure(Ωc,2)
@test sum(∫(Operation(det)(g)/sqrt_abs_g)dΩc) ≈ 1.5*3.0 - 0.09

# A triangulation of part of the model reads the metric of the cells it holds
Ωp = Triangulation(cmodel,[2,4])
gp = first_point_values(first(get_metric_data(Ωp)))
@test [v[1,1] for v in gp] == expected_g_11[[2,4]]

# A facet needs the induced metric, which is a different object
@test_throws ErrorException get_metric_data(BoundaryTriangulation(cmodel))
@test_throws ErrorException get_metric_data(SkeletonTriangulation(cmodel))

# A model whose cells are embedded in a space of higher dimension is rejected:
# its physical points carry the dimension of that space, so its cell fields
# cannot be a metric of the cells' own coordinates.
flat_quad = UnstructuredGrid([Point(0.,0.,0.),Point(1.,0.,0.),Point(0.,1.,0.),Point(1.,1.,0.)],
                             Table([[1,2,3,4]]), [LagrangianRefFE(Float64,QUAD,1)], [1])
embedded = UnstructuredDiscreteModel(flat_quad)
@test (num_cell_dims(embedded), num_point_dims(embedded)) == (2,3)
@test get_metric_data(Triangulation(embedded)) === (LinearAlgebra.I,)
@test_throws ErrorException get_metric_data(Triangulation(chart_model(embedded,CovChartMetric)))

# ── The variance of the metric is the model's to carry ────────────────────────

# An un-varianced metric is rejected rather than read as covariant
@test_throws ArgumentError get_metric_data(Triangulation(chart_model(model,ChartMetric)))

# A contravariant tensor is an inverse metric, not a metric
@test_throws ArgumentError get_metric_data(Triangulation(chart_model(model,ConChartMetric)))

# The same holds of an overridden g⁻¹, which must be contravariant
unvarianced_inv = ChartModel(model,
                             cell_metric_fields(CovChartMetric,model),
                             cell_metric_fields(ChartMetric,model))
@test_throws ArgumentError get_metric_data(Triangulation(unvarianced_inv))

# A model holding a closed form for g⁻¹ keeps it rather than inverting g
imodel = ChartModel(model,
                    cell_metric_fields(CovChartMetric,model),
                    cell_metric_fields(InvChartMetric,model))
@test get_cell_inv_metric(imodel) === imodel.cell_inv_metric
g, g_inv = get_metric_data(Triangulation(imodel))
gv, giv = map(first_point_values, (g,g_inv))
@test all(v -> v isa SymTensorValue{Con(2),Float64}, giv)
@test all(cell -> get_array(giv[cell]) ≈ inv(get_array(gv[cell])), 1:num_cells(imodel))

end # module
