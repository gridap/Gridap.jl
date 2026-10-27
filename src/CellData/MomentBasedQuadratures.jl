# Moment-based quadratures and measures on triangulations.

# Patch gathering

"""
    PatchMeasure(ctrian::Triangulation, fmeas::Measure, ptopo::PatchTopology)
    PatchMeasure(ctrian::Triangulation, fmeas::Measure, ctface_to_ftfaces)

Returns a measure on `ctrian` obtained by gathering, for each cell of `ctrian`, the
quadrature points of `fmeas` on the fine cells it contains. See [`PatchQuadrature`](@ref).
"""
function PatchMeasure(ctrian::Triangulation, fmeas::Measure, patches)
  quads = PatchQuadrature(ctrian, get_cell_quadrature(fmeas), patches)
  Measure(quads)
end

"""
    PatchQuadrature(ctrian::Triangulation, fquads::CellQuadrature, ptopo::PatchTopology)

Same as `PatchQuadrature(ctrian, fquads, ctface_to_ftfaces)`, with the table built from `ptopo`.
The cells of the background model of `ctrian` are the patches of `ptopo`, and the
cells of the background model of `get_triangulation(fquads)` are the cells of `ptopo`.
Fine cells that are not in the fine triangulation are skipped.
"""
function PatchQuadrature(ctrian::Triangulation, fquads::CellQuadrature, ptopo::PatchTopology)
  Df = num_cell_dims(ctrian)
  ftrian = get_triangulation(fquads)
  patch_to_fcells = Geometry.get_patch_cells(ptopo)
  fcell_to_ftface = Geometry.get_glue(ftrian,Val(Df)).mface_to_tface
  if isa(fcell_to_ftface, IdentityVector)
    patch_to_ftfaces = patch_to_fcells
  else
    pface_to_ftface = collect(Int32, Geometry.patch_reindex(ptopo, fcell_to_ftface))
    patch_to_ftfaces = Table(map(
      ftfaces -> filter(>(0), ftfaces), Table(pface_to_ftface, patch_to_fcells.ptrs)
    ))
  end
  ctface_to_patch = Geometry.get_glue(ctrian,Val(Df)).tface_to_mface
  ctface_to_ftfaces = lazy_map(Reindex(patch_to_ftfaces), ctface_to_patch)
  PatchQuadrature(ctrian, fquads, ctface_to_ftfaces)
end

"""
    PatchQuadrature(ctrian::Triangulation, fquads::CellQuadrature, ctface_to_ftfaces)

Returns a `CellQuadrature` on `ctrian` where each cell `c` carries the quadrature points
and weights of `fquads` on the cells `ctface_to_ftfaces[c]` of `get_triangulation(fquads)`,
in the physical domain. Both cells are indexed by triangulation (not model) ids.
"""
function PatchQuadrature(ctrian::Triangulation, fquads::CellQuadrature, ctface_to_ftfaces::AbstractArray)
  Df = num_cell_dims(ctrian)
  ftrian = get_triangulation(fquads)
  @check num_cell_dims(ftrian) == Df
  @check length(ctface_to_ftfaces) == num_cells(ctrian)

  # Move to the physical domain
  if DomainStyle(fquads) == PhysicalDomain()
    ftface_points = fquads.cell_point
    ftface_weights = fquads.cell_weight
  else
    ftface_points = get_data(change_domain(get_cell_points(fquads), PhysicalDomain()))
    ftface_Jt = lazy_map(evaluate, lazy_map(∇, get_cell_map(ftrian)), fquads.cell_point)
    ftface_weights = lazy_map(Broadcasting((w,Jt) -> w*meas(Jt)), fquads.cell_weight, ftface_Jt)
  end

  # From the fine triangulation to the coarse triangulation
  ctface_points = lazy_map(Arrays.AppendEntriesMap(ftface_points), ctface_to_ftfaces)
  ctface_weights = lazy_map(Arrays.AppendEntriesMap(ftface_weights), ctface_to_ftfaces)
  ctface_quads = lazy_map(GenericQuadrature, ctface_points, ctface_weights, Fill("Patch Quadrature", num_cells(ctrian)))

  return CellQuadrature(ctface_quads, ctface_points, ctface_weights, ctrian, PhysicalDomain())
end

# Moment-based measures

"""
    MomentBasedMeasure(trian::Triangulation, meas::Measure, degree; T=Float64, space=:Q)
    MomentBasedMeasure(trian::Triangulation, fmeas::Measure, patches, degree; kwargs...)

Returns a measure on `trian` where each cell carries a moment-based quadrature (see
[`MomentBasedQuadratureMap`](@ref)), exact for the polynomial space `space ∈ (:P,:Q)` of
order `degree`. Two constructors are provided:

- If a `meas::Measure` is given, that can integrate polynomials on each cell of `trian`,
  the moment-based quadrature is constructed directly from it.

- Alternatively, one can provide a measure `fmeas` on a fine triangulation along with
  `patches`, providing how the cells of `trian` are built from the fine cells of `fmeas`.
  `patches` is either a `PatchTopology` whose patches are the cells of the background model
  of `trian` (agglomeration), or a table `ctface_to_ftfaces` giving, for each cell of `trian`,
  the fine cells of `fmeas` it is built from (e.g. the sub-cells of a cut cell).
  The moments are then computed by gathering the quadratures of `fmeas` on the fine cells
  within each cell of `trian` (see [`PatchQuadrature`](@ref)).
"""
function MomentBasedMeasure(trian::Triangulation, meas::Measure, degree; kwargs...)
  quad = MomentBasedQuadrature(trian, get_cell_quadrature(meas), degree; kwargs...)
  Measure(quad)
end

function MomentBasedMeasure(trian::Triangulation, fmeas::Measure, patches, degree; kwargs...)
  cmeas = PatchMeasure(trian, fmeas, patches)
  MomentBasedMeasure(trian, cmeas, degree; kwargs...)
end

"""
    MomentBasedQuadrature(trian::Triangulation, quad::CellQuadrature, degree; T=Float64, space=:Q)

Returns a `CellQuadrature` on `trian` where each cell carries a moment-based quadrature,
with moments computed with `quad`. The result lives in the same domain as `quad`:
if `quad` is in the physical domain, the polynomial space is defined on the physical
bounding box of each cell. Otherwise, it is defined on the bounding box of the cell polytope.
"""
function MomentBasedQuadrature(
  trian::Triangulation{Dc}, quad::CellQuadrature, degree;
  T::Type = Float64, space::Symbol = :Q
) where Dc
  k = MomentBasedQuadratureMap(Dc, degree; T=T, space=space)
  data_domain_style = DomainStyle(quad)
  if data_domain_style == PhysicalDomain()
    cell_boxes = get_cell_coordinates(trian)
  else
    cell_boxes = Geometry.get_cell_polytopes(trian)
  end
  cell_quads = Arrays.lazy_collect(lazy_map(k, cell_boxes, get_data(quad)))
  return CellQuadrature(trian, cell_quads, data_domain_style, PhysicalDomain())
end
