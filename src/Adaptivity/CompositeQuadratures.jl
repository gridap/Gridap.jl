
function CellData.CellQuadrature(
  trian::Triangulation,
  rrules::AbstractVector{<:RefinementRule},
  degree::Integer;
  kwargs...
)
  return CellData.CellQuadrature(trian,rrules,degree,PhysicalDomain();kwargs...)
end

function CellData.CellQuadrature(
  trian::Triangulation,
  rrules::AbstractVector{<:RefinementRule},
  degree::Integer,
  ids::DomainStyle;
  kwargs...
)
  cell_quad = lazy_map(rr -> Quadrature(rr,degree;kwargs...),rrules)
  return CellData.CellQuadrature(trian,cell_quad,integration_domain_style=ids)
end

struct CompositeQuadrature <: QuadratureName end

function ReferenceFEs.Quadrature(rr::RefinementRule,degree::Integer;kwargs...)
  return ReferenceFEs.Quadrature(ReferenceFEs.get_polytope(rr),CompositeQuadrature(),rr,degree;kwargs...)
end

function ReferenceFEs.Quadrature(
  p::Polytope,::CompositeQuadrature,rr::RefinementRule,degree::Integer;kwargs...
)
  @check p === ReferenceFEs.get_polytope(rr)
  subgrid  = get_ref_grid(rr)
  cellmap  = get_cell_map(rr)
  measures = get_cell_measures(rr)
  quads    = ReferenceFEs.Quadrature(subgrid,degree;kwargs...)

  npts = sum(map(num_points,quads))
  WT = eltype(get_weights(first(quads)))
  CT = eltype(get_coordinates(first(quads)))
  weights = Vector{WT}(undef,npts)
  cpoints = Vector{CT}(undef,npts)
  conn = Vector{Vector{Int32}}(undef,length(quads))
  k = 1
  for (iq,q) in enumerate(quads)
    n = num_points(q)
    w = get_weights(q)
    c = get_coordinates(q)
    weights[k:k+n-1] .= w .* measures[iq]
    cpoints[k:k+n-1] .= (cellmap[iq])(c)
    conn[iq] = collect(k:k+n-1)
    k = k+n
  end
  fpoints = map(get_coordinates,quads)
  ids = FineToCoarseIndices(conn)
  coordinates = FineToCoarseArray{eltype(cpoints)}(rr,cpoints,fpoints,ids)
  return GenericQuadrature(coordinates,weights,"Composite quadrature")
end
