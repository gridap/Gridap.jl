# The trace of a purely covariant tensor field, a differential form in
# particular, is the pullback of the inclusion, and nothing else. Only a purely
# covariant tensor has one, see [`_check_purely_covariant`](@ref).
#
# With `D` the dimension of the background `Ω` and `Dm = D-1`
# that of `Γ`, the maps at play are
#
#   φk  : K̂ → K ⊂ ℝᴰ   reference cell to physical cell,    `Jk  = ∂φk`  (D ×D)
#   ιrf : F̂ → K̂        reference face in reference cell,   `Jrf = ∂ιrf` (Dm×D)
#   φf  : F̂ → F ⊂ ℝᴰ   reference face to physical face,    `Jf  = ∂φf`  (Dm×D)
#
# with `φf = φk ∘ ιrf`, hence `Jf = Jrf ⋅ Jk`.
#
# `ι : Γ → Ω` is the identity on physical points, so its differential is the
# inclusion of `TₓΓ` into `TₓΩ` and has no components until a frame is chosen
# on each side. Read in the ambient frame of `Ω` and the reference face frame
# of `Γ`, which are the two Gridap evaluates fields in, it is exactly `Jf`:
# there is no inclusion Jacobian separate from the face map one, and
# `pullback(ω, Jf)` is already `ι*ω`, with no reference to a normal vector.
#
# `Jf` is a frame, not a metric and not an orientation, so it carries neither.
# The metric `Γ` inherits is handled by [`get_metric_data`](@ref) and the
# orientation `Ω` induces on it by the integration of its top forms.

struct TensorTrace{A<:CellField}
  dι::A
end

# applied like the lifted operators of a weak form, the differential reaching
# the quadrature points as an operand of `pullback`. The guard is here and not
# on the constructor, which only sees `dι`.
function (tr::TensorTrace)(a::CellField)
  _check_purely_covariant(_cell_field_value_type(a))
  pullback(a, tr.dι)
end

# `dι` alone carries a covariant index, so `ι*` restricts a purely covariant
# tensor and only one: a contravariant index would need a pseudo inverse of
# `dι`, which projects the ambient frame onto `TΓ` rather than restricting
# anything, and an un-varianced one has no transformation rule at all.
_check_purely_covariant(::Type) = nothing

function _check_purely_covariant(::Type{T}) where T<:MultiValue
  all(D -> D isa Cov, size(T)) || throw(ArgumentError("""\n
  The trace of a tensor field is the pullback of the inclusion, which is defined
  for a purely covariant tensor, of type `(0,s,0)`. Got `$T`, whose indices
  `$(size(T))` are not all covariant.
  """))
  nothing
end

# the value type a cell field takes, through the arrays a basis evaluates to
_cell_field_value_type(a::CellField) =
  _integrand_value_type(typeof(testitem(a(get_cell_points(get_triangulation(a))))))

"""
    _get_inclusion_jac(Γ::Triangulation)

`dι` read in the ambient frame of the background and the reference cell frame
of `Γ`, that is the transposed Jacobian of the cell map of `Γ`, as the
`CellField` of `Dm×D` tensors that [`pullback`](@ref) takes as a second
argument. It is the Jacobian of an inclusion only for a `Γ` embedded in its
background, which the public [`get_trace_operator`](@ref) is what checks: hence
internal.
"""
_get_inclusion_jac(Γ::Triangulation) =
  GenericCellField(lazy_map(Broadcasting(∂), get_cell_map(Γ)), Γ, ReferenceDomain())

"""
    get_trace_operator(Γ::Triangulation, Ω::Triangulation)

The trace operator of a lower dimensional triangulation `Γ` of the background
of `Ω`: `trΓ(ω)` is the pullback `ι*ω` along the inclusion `ι : Γ → Ω`, what a
tensor field of `Ω` restricts to on `Γ`. Only a purely covariant field, of type
`(0,s,0)`, has such a restriction, a differential form in particular.

    trΓ = get_trace_operator(Γ, Ω)
    ∫( p̄ ∧ trΓ(v) )dΓ       # or, as a function, ∫( p̄ ∧ ⋆(trΓ(v)) )dΓ

The traced field has its components in the reference cell coordinates of `Γ`,
where the metric is [`get_metric_data`](@ref)`(Γ)` and a top form integrates
with the orientation the background induces, so that Stokes' theorem

    ∫_Ω 𝑑ω = ∮_{∂Ω} trΓ(ω)

holds in every degree -- unlike a trace defined by a rescaled Jacobian, which
can only be made exact for one degree at a time.
"""
function get_trace_operator(Γ::Triangulation, Ω::Triangulation)
  @check Geometry.is_change_possible(Ω,Γ) """\n
  `Γ` must be a triangulation of faces of `Ω` to be included into it.
  """
  @check num_cell_dims(Γ) < num_cell_dims(Ω) """\n
  The trace of a form is taken on a triangulation of lower dimension than the
  background it is a form of, got $(num_cell_dims(Γ)) and $(num_cell_dims(Ω)).
  """
  TensorTrace(_get_inclusion_jac(Γ))
end

# ── The orientation the background induces ────────────────────────────────────

# A facet is parametrised by its own reference map, ordered by its vertex
# numbering, which has nothing to do with the orientation its background
# induces, the one making `det(n, t₁, …, t_{D-1})` positive for `n` the outward
# normal. The two differ by a sign per facet,
#
#   ε = sign(det Jk) ⋅ ε_ref(nref, Jrf),
#
# `ε_ref` comparing the facet's parametrisation to the induced one on the
# reference cell, and `sign(det Jk)` carrying the cell's own orientation over to
# physical space, which `simplexify` makes non constant. Both are reference or
# cell-wise data: no normal field and no quadrature point is involved, and the
# sign cannot be folded into the Jacobian, where a scalar would act on a degree
# `p` pullback as its `p`-th power.
#
# `Jrf` is the differential of `ιrf`, the map the glue carries from the
# reference facet into the reference cell, which is the facet's parametrisation
# itself: reading the tangents off it assumes nothing of how a facet is
# parametrised and takes its vertex permutation into account on its own.

"""
    _reference_induced_sign(nref::VectorValue, Jrf::MultiValue)

`+1` if the parametrisation of a facet of transposed Jacobian `Jrf` into its
reference cell is positively oriented for the orientation that cell induces on
it, `-1` otherwise, `nref` being the outward normal of the facet in the
reference cell.
"""
function _reference_induced_sign(
  nref::VectorValue{D}, Jrf::MultiValue{Tuple{Dm,Dn}}) where {D,Dm,Dn}

  # rows of `M`: the outward normal, then the tangents of the facet's map
  M = TensorValue{D,D}(ntuple(Val(D*D)) do k
    i, j = mod1(k,D), cld(k,D)
    i == 1 ? nref[j] : Jrf[i-1,j]
  end)
  sign(det(M))
end

# `sign(det Jk)`, one value per cell of the background. Taken at one point of
# each cell: a cell map whose Jacobian changes sign inside a cell is degenerate.
function _cell_orientation(model::DiscreteModel)
  Ω  = Triangulation(model)
  Jk = GenericCellField(lazy_map(Broadcasting(∂), get_cell_map(Ω)), Ω, ReferenceDomain())
  a  = Jk(get_cell_points(Ω))
  c  = array_cache(a)
  [sign(det(first(getindex!(c,a,cell)))) for cell in 1:num_cells(Ω)]
end

"""
    induced_orientation(Γ::BoundaryTriangulation)

The sign, one per facet of `Γ`, relating the facet's own reference
parametrisation to the orientation the background induces on it.
"""
function induced_orientation(Γ::BoundaryTriangulation)
  glue    = Γ.glue
  model   = get_background_model(Γ)
  cell_or = _cell_orientation(model)
  ctype_lface_nref = [get_facet_normal(p) for p in get_polytopes(get_grid(model))]
  # `ιrf` is affine in the facet coordinates, so its differential is constant
  # over a facet and one point of each tells its orientation.
  Jrf = GenericCellField(
    lazy_map(Broadcasting(∂), get_glue(Γ, Val(num_cell_dims(model))).tface_to_mface_map),
    Γ, ReferenceDomain())
  a = Jrf(get_cell_points(Γ))
  c = array_cache(a)
  map(1:num_cells(Γ)) do face
    cell = glue.face_to_cell[face]
    nref = ctype_lface_nref[glue.cell_to_ctype[cell]][glue.face_to_lface[face]]
    cell_or[cell] * _reference_induced_sign(nref, first(getindex!(c,a,face)))
  end
end

# The orientation the cells of a triangulation carry, one sign per cell, as
# [`get_metric_data`](@ref) and the integral of a top form need it. Unlike the
# metric it is not a function of the background, which is why it is the one
# thing a new triangulation has to declare.
_induced_orientation(trian::Triangulation) = @abstractmethod """\n
The orientation the cells of a `$(nameof(typeof(trian)))` carry is not
implemented, and there is no default: a triangulation that is a boundary of its
background takes the orientation that background induces on it, one that is not
takes its own parametrisation, an interior facet takes the one its `plus` side
induces, and a wrong sign between any two of those is silent.
"""

# nothing induces an orientation on these, their cells are oriented as they are
# parametrised
_induced_orientation(trian::BodyFittedTriangulation) = Fill(1.0, num_cells(trian))
_induced_orientation(trian::GenericTriangulation)    = Fill(1.0, num_cells(trian))

_induced_orientation(Γ::BoundaryTriangulation) = induced_orientation(Γ)

# The two cells an interior facet separates induce opposite orientations on it,
# so one side has to name it. `Λ` is parametrised by its `plus` side,
# `get_grid(Λ) = get_grid(Λ.plus)`, and `get_normal_vector(Λ).⁺` points out of
# that side, so its orientation is `[ι_N μ]` for that `N` -- the choice `jump`
# makes already. With `jump(a) = a.⁺ - a.⁻` each cell `K` then sees the facet
# with the orientation it induces, which is what a cell-wise Stokes' theorem
# asks for:
#
#   Σ_K ∮_{∂K} ι*ω = ∫_Λ jump(ι*ω) + ∮_{∂Ω} ι*ω.
_induced_orientation(Λ::SkeletonTriangulation) = _induced_orientation(Λ.plus)

# the cells of these are those of another triangulation, and so is their
# orientation
_induced_orientation(trian::Geometry.TriangulationView) =
  lazy_map(Reindex(_induced_orientation(trian.parent)), trian.cell_to_parent_cell)

_induced_orientation(trian::AppendedTriangulation) =
  lazy_append(_induced_orientation(trian.a), _induced_orientation(trian.b))

# a composite is its `dtrian` placed on its `rtrian`, whose cells it has
_induced_orientation(trian::Geometry.CompositeTriangulation) =
  _induced_orientation(trian.dtrian)

# a patch triangulation relabels the patches, not the cells: `get_grid` is the
# grid of `trian`, so a patch boundary carries the orientation of the boundary
# it is built of
_induced_orientation(trian::PatchTriangulation) = _induced_orientation(trian.trian)
