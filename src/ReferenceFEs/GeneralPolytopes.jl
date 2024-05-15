"""
    struct GeneralPolytope{D,Dp,Tp} <: Polytope{D}

  The `GeneralPolytope` is definded defined by a set of vertices and a rototation
  system (a planar oriented graph). This polytopal representation can represent
  any polytope in 2 and 3 dimensions.

  In 2 dimensions ([`Polygon`](@ref)), the representation of the polygon is a closed polyline.

  In 3 dimensions ([`Polyhedron`](@ref)), the rotation system generates the connectivities, each   facet is a closed cycle of the graph.
  This construction allows complex geometrical operations, e.g., intersecting polytopes by halfspaces.
  See also,

  > K. Sugihara, "A robust and consistent algorithm for intersecting convex polyhedra", Comput. Graph. Forum 13 (3) (1994) 45–54, doi: [10.1111/1467-8659.1330045](https://doi.org/10.1111/1467-8659.1330045)

  > D. Powell, T. Abel, "An exact general remeshing scheme applied to physically conservative voxelization", J. Comput. Phys. 297 (Sept. 2015) 340–356, doi: [10.1016/j.jcp.2015.05.022](https://doi.org/10.1016/j.jcp.2015.05.022.

  > S. Badia, P. A. Martorell, F. Verdugo. "Geometrical discretisations for unfitted finite elements on explicit boundary representations", J.Comput. Phys. 460 (2022): 111162. doi: [10.1016/j.jcp.2022.111162](https://doi.org/10.1016/j.jcp.2022.111162)
"""
struct GeneralPolytope{D,Dp,Tp,Td} <: Polytope{D}
  vertices::Vector{Point{Dp,Tp}}
  edge_vertex_graph::Vector{Vector{Int32}}
  n_m_to_nface_to_mfaces::Matrix{Vector{Vector{Int}}}
  facedims::Vector{Int32}
  dimranges::Vector{UnitRange{Int}}
  dface_nfaces::Vector{Vector{Int}}
  isopen::Bool
  metadata::Td
  function GeneralPolytope{D}(
    vertices::Vector{Point{Dp,Tp}},
    edge_vertex_graph::Vector{Vector{Int32}},
    n_m_to_nface_to_mfaces::Matrix{Vector{Vector{Int}}},
    facedims::Vector{Int32},
    dimranges::Vector{UnitRange{Int}},
    dface_nfaces::Vector{Vector{Int}},
    isopen::Bool,
    metadata::Td) where {D,Dp,Tp,Td}
    new{D,Dp,Tp,Td}(
      vertices,
      edge_vertex_graph,
      n_m_to_nface_to_mfaces,
      facedims,
      dimranges,
      dface_nfaces,
      isopen,
      metadata)
  end
end

"""
    Polygon = GeneralPolytope{2}

  A polygon is a [`GeneralPolytope`](@ref) in 2 dimensions.
"""
const Polygon = GeneralPolytope{2}

"""
    Polyhedron = GeneralPolytope{3}

  A polyhedron is a [`GeneralPolytope`](@ref) in 3 dimensions.
"""
const Polyhedron = GeneralPolytope{3}

# Constructors

function GeneralPolytope{D}(
  vertices::Vector{<:Point},
  graph::Vector{Vector{Int32}},
  isopen::Bool,
  data) where D

  n = D+1
  n_m_to_nface_to_mfaces = Matrix{Vector{Vector{Int}}}(undef,n,n )
  dimranges = Vector{UnitRange{Int}}(undef,0)
  dface_nfaces = Vector{Vector{Int}}(undef,0)
  facedims = Vector{Int32}(undef,0)
  GeneralPolytope{D}(
    vertices,
    graph,
    n_m_to_nface_to_mfaces,
    facedims,
    dimranges,
    dface_nfaces,
    isopen,
    data)
end

function GeneralPolytope{D}(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector},
  isopen::Bool,
  data) where D

  GeneralPolytope{D}(collect(vertices),graph,isopen,data)
end

"""
    GeneralPolytope{D}(vertices,graph;kwargs...)

  Constructor of a [`GeneralPolytope`](@ref) that generates a polytope of
  D dimensions with the given `vertices` and `graph` of connectivities.
"""
function GeneralPolytope{D}(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector}
  ;isopen=false::Bool,
  metadata=nothing) where D

  GeneralPolytope{D}(vertices,graph,isopen,metadata)
end

function GeneralPolytope{D}(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector},
  data) where D

  isopen = false
  GeneralPolytope{D}(vertices,graph,isopen,data)
end

function generate_polytope_data(p::Polytope;metadata=nothing)
  generate_polytope_data(p,metadata)
end

function generate_polytope_data(p::Polytope,::Nothing)
  nothing
end

function Polygon(p::Polytope{2},vertices::AbstractVector{<:Point};kwargs...)
  if p == TRI
    e_v_graph = [[2,3],[3,1],[1,2]]
    perm = [1,2,3]
  elseif p == QUAD
    e_v_graph = [[2, 3],[4, 1],[1, 4],[3, 2]]
    perm = [1,2,4,3]
  else
    @unreachable
  end
  vertices = map(Reindex(vertices),perm)
  e_v_graph = map(Reindex(e_v_graph),perm)
  e_v_graph = map(i->replace(i, Dict(perm .=> 1:length(perm))...),e_v_graph)
  e_v_graph = map(i->Int32.(i),e_v_graph)
  data = generate_polytope_data(p;kwargs...)
  Polygon(vertices,e_v_graph,data)
end

function Polygon(vertices::AbstractVector{<:Point};kwargs...)
  graph = map(1:length(vertices)) do i
    inext = i == length(vertices) ? 1 : i+1
    iprev = i == 1 ? length(vertices) : i-1
    Int32[iprev,inext]
  end
  Polygon(vertices,graph;kwargs...)
end

function Polyhedron(p::Polytope{3},vertices::AbstractVector{<:Point};kwargs...)
  if p == TET
    e_v_graph = [[2,4,3],[3,4,1],[1,4,2],[1,2,3]]
  elseif p == HEX
    e_v_graph = [
      [5, 2, 3],
      [6, 4, 1],
      [7, 1, 4],
      [8, 3, 2],
      [1, 7, 6],
      [2, 5, 8],
      [3, 8, 5],
      [4, 6, 7] ]
  else
    @unreachable
  end
  e_v_graph = map(i->Int32.(i),e_v_graph)
  data = generate_polytope_data(p;kwargs...)
  Polyhedron(vertices,e_v_graph,data)
end

function GeneralPolytope{D}(p::Polytope;kwargs...) where D
  GeneralPolytope{D}(p,get_vertex_coordinates(p);kwargs...)
end

# Interface

num_dims(::T) where T<:GeneralPolytope = num_dims(T)

num_dims(::Type{<:GeneralPolytope{D}}) where D = D

num_cell_dims(a::GeneralPolytope) = num_dims(a)

point_eltype(::T) where T<:GeneralPolytope = point_eltype(T)

point_eltype(::Type{<:GeneralPolytope{D,Dp,T}}) where {D,Dp,T} = T

num_point_dims(::Type{<:GeneralPolytope{D,Dp}}) where {D,Dp} = Dp

num_vertices(a::GeneralPolytope) = length(a.vertices)

get_vertex_coordinates(a::GeneralPolytope) = a.vertices

Base.getindex(a::GeneralPolytope,i::Integer) = a.vertices[i]

"""
    get_graph(p::GeneralPolytope) -> Vector{Vector{Int32}}

  It returns the edge-vertex graph of the polytope `p`.
"""
@inline get_graph(a::GeneralPolytope) = a.edge_vertex_graph


"""
    get_metadata(p::GeneralPolytope)

  It return the metadata stored in the polytope `p`.
"""
get_metadata(a::GeneralPolytope) = a.metadata


"""
    isopen(p::GeneralPolytope) -> Bool

  In return whether the polytope is watter tight or not.
"""
Base.isopen(a::GeneralPolytope) = a.isopen

"""
    isactive(p::GeneralPolytope,vertex::Integer) -> Bool

  It returns whether a vertex is connected to any other vertex.
"""
function isactive(p::Polyhedron,vertex::Integer)
  !isempty( get_graph(p)[vertex] )
end

"""
    check_polytope_graph(p::GeneralPolytope) -> Bool

  It checks whether the graph is well-constructed. The graph must be oriented
  and planar.
"""
function check_polytope_graph(p::GeneralPolytope)
  check_polytope_graph(get_graph(p))
end

function check_polytope_graph(graph::AbstractVector{<:AbstractVector})
  for v in 1:length(graph)
    !isempty(graph[v]) || continue
    for vneig in graph[v]
      vneig > 0 || continue
      v ∈ graph[vneig] || return false
    end
  end
  true
end

is_simplex(::GeneralPolytope) = false

is_n_cube(::GeneralPolytope) = false

function simplexify(p::GeneralPolytope{D}) where D
  @assert !isopen(p)
  X,T = simplexify_interior(p)
  @check X == get_vertex_coordinates(p)
  T, simplex_polytope(Val{D}())
end

simplex_polytope(::Val{0}) = VERTEX

simplex_polytope(::Val{1}) = SEGMENT

simplex_polytope(::Val{2}) = TRI

simplex_polytope(::Val{3}) = TET

function Polytope{0}(p::GeneralPolytope,faceid::Integer)
  VERTEX
end

function Polytope{1}(p::GeneralPolytope,faceid::Integer)
  SEGMENT
end

function Polytope{D}(p::GeneralPolytope{D},faceid::Integer) where D
  p
end

function Polytope{2}(p::Polyhedron,faceid::Integer)
  f_to_v = get_faces(p,2,0)
  coords = get_vertex_coordinates(p)
  vertices = coords[f_to_v[faceid]]
  Polygon(vertices)
end

Base.:(==)(a::GeneralPolytope,b::GeneralPolytope) = false

function Base.:(==)(a::GeneralPolytope{D},b::GeneralPolytope{D}) where D
  a === b ||
  ( get_vertex_coordinates(a) == get_vertex_coordinates(b) &&
  get_graph(a) == get_graph(b) )
end

function num_faces(p::GeneralPolytope{D},d::Integer) where D
  if d == 0
    num_vertices(p)
  elseif d == D
    1
  else
    length( get_faces(p,d,0) )
  end
end

function get_faces(p::GeneralPolytope,n::Integer,m::Integer)
  setup_faces!(p,n,m)
  p.n_m_to_nface_to_mfaces[n+1,m+1]
end

function get_facet_orientations(p::GeneralPolytope)
  ones(Int,num_facets(p))
end

function get_facet_normal(p::Polyhedron)
  D = 3
  f_to_v = get_faces(p,D-1,0)
  coords = get_vertex_coordinates(p)
  map(f_to_v) do v
    v1 = coords[v[2]]-coords[v[1]]
    v2 = coords[v[3]]-coords[v[1]]
    v1 /= norm(v1)
    v2 /= norm(v2)
    n = v1 × v2
    n /= norm(n)
  end
end

function get_facet_normal(p::Polygon)
  D = 2
  f_to_v = get_faces(p,D-1,0)
  coords = get_vertex_coordinates(p)
  @notimplementedif num_components(eltype(coords)) != 2
  map(f_to_v) do v
    e = coords[v[2]]-coords[v[1]]
    n = VectorValue( e[2], -e[1] )
    n /= norm(n)
  end
end

function get_edge_tangent(p::GeneralPolytope)
  e_to_v = get_faces(p,1,0)
  coords = get_vertex_coordinates(p)
  map(e_to_v) do v
    e = coords[v[2]]-coords[v[1]]
    e / norm(e)
  end
end

function get_dimranges(p::GeneralPolytope)
  setup_dimranges!(p)
  p.dimranges
end

function get_dimrange(p::GeneralPolytope,d::Integer)
  setup_dimranges!(p)
  p.dimranges[d+1]
end

function get_faces(p::GeneralPolytope)
  setup_face_to_nfaces!(p)
  p.dface_nfaces
end

function get_facedims(p::GeneralPolytope)
  setup_facedims!(p)
  p.facedims
end

function setup_dimranges!(p::GeneralPolytope{D}) where D
  if length(p.dimranges) < D+1
    lens = map(i->num_faces(p,i),0:D)
    sums = cumsum(lens)
    resize!(p.dimranges,D+1)
    for (i,(l,s)) in enumerate(zip(lens,sums))
      p.dimranges[i] = s-l+1 : s
    end
  end
end

function setup_face_to_nfaces!(p::GeneralPolytope)
  if length(p.dface_nfaces) < num_vertices(p)
    facedims = get_facedims(p)
    dface_nfaces = Vector{Vector{Int}}(undef,length(facedims))
    ofsets = get_offsets(p)
    for (f,d) in enumerate(facedims)
      df = f - ofsets[d+1]
      nfs = Int[]
      for n in 0:d
        union!(nfs,get_faces(p,d,n)[df] .+ ofsets[n+1])
      end
      dface_nfaces[f] = nfs
    end
    copy!(p.dface_nfaces,dface_nfaces)
  end
  nothing
end

function setup_facedims!(p::GeneralPolytope)
  if length(p.facedims) < num_vertices(p)
    dimranges = get_dimranges(p)
    n_faces = dimranges[end][end]
    facedims = _nface_to_nfacedim(n_faces,dimranges)
    copy!(p.facedims,facedims)
  end
end

function setup_faces!(p::GeneralPolytope{D},dimfrom,dimto) where D
  if isassigned(p.n_m_to_nface_to_mfaces,dimfrom+1,dimto+1)
    return nothing
  end
  if dimfrom == dimto
    setup_nface_to_nface!(p,dimfrom)
  elseif dimfrom == D
    setup_cell_to_faces!(p,dimto)
  elseif dimto == 0
    setup_face_to_vertices!(p,dimfrom)
  elseif dimfrom > dimto
    setup_nface_to_mface!(p,dimfrom,dimto)
  else
    setup_nface_to_mface_dual!(p,dimto,dimfrom)
  end
  nothing
end

function setup_face_to_vertices!(p::GeneralPolytope,d)
  if !isassigned(p.n_m_to_nface_to_mfaces,d+1,1)
    df_to_v = generate_face_to_vertices(p,d)
    p.n_m_to_nface_to_mfaces[d+1,1] = df_to_v
  end
end

function setup_cell_to_faces!(p::GeneralPolytope{D},d) where D
  if !isassigned(p.n_m_to_nface_to_mfaces,D+1,d+1)
    num_f = num_faces(p,d)
    c_to_f = [ collect(1:num_f) ]
    p.n_m_to_nface_to_mfaces[D+1,d+1] = c_to_f
  end
end

function setup_nface_to_nface!(p::GeneralPolytope,n)
  if !isassigned(p.n_m_to_nface_to_mfaces,n+1,n+1)
    num_nf = num_faces(p,n)
    nf_to_nf = map(i->Int[i],1:num_nf)
    p.n_m_to_nface_to_mfaces[n+1,n+1] = nf_to_nf
  end
end

function setup_nface_to_mface!(p::Polyhedron,n,m)
  if !isassigned(p.n_m_to_nface_to_mfaces,n+1,m+1)
    @notimplementedif n != 2 && m != 1
    nf_to_v = get_faces(p,n,0)
    v_to_mf = get_faces(p,0,m)
    nf_to_ftype = map( length, nf_to_v )
    ftype_to_lmf_to_lv = map(1:maximum(nf_to_ftype)) do ftype
      map(1:ftype) do i
        inext = i == ftype ? 1 : i+1
        Int[i,inext]
      end
    end
    nf_to_mf = Gridap.Geometry.find_cell_to_faces(
      Table(nf_to_v),
      ftype_to_lmf_to_lv,
      nf_to_ftype,
      Table(v_to_mf))
    p.n_m_to_nface_to_mfaces[n+1,m+1] = Vector(nf_to_mf)
  end
end

function setup_nface_to_mface_dual!(p::GeneralPolytope,dimto,dimfrom)
  if !isassigned(p.n_m_to_nface_to_mfaces,dimfrom+1,dimto+1)
    @assert dimfrom < dimto
    nf_to_mf = get_faces(p,dimto,dimfrom)
    n_mf = num_faces(p,dimfrom)
    mf_to_nf = Gridap.Geometry.generate_cells_around(Table(nf_to_mf),n_mf)
    p.n_m_to_nface_to_mfaces[dimfrom+1,dimto+1] = Vector(mf_to_nf)
  end
end

function generate_face_to_vertices(p::Polyhedron,d::Integer)
  if d == 1
    generate_edge_to_vertices(p)
  elseif d == 2
    generate_facet_to_vertices(p)
  else
    @unreachable
  end
end

function generate_face_to_vertices(p::Polygon,d::Integer)
  if d == 1
    generate_facet_to_vertices(p)
  else
    @unreachable
  end
end

function generate_facet_to_vertices(poly::Polyhedron)
  D = 3
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  T = Vector{Int32}[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      vnext > 0 || continue
      k = [v]
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        vnext > 0 || break
        push!(k,vcurrent)
      end
      if length(k) >=D
        push!(T,k)
      end
    end
  end
  T
end

function generate_edge_to_vertices(poly::GeneralPolytope)
  graph = get_graph(poly)
  T = Vector{Int32}[]
  for v in 1:length(graph)
    for vneig in graph[v]
      if vneig > v
        push!(T,[v,vneig])
      end
    end
  end
  T
end

function generate_facet_to_vertices(poly::Polygon)
  graph = get_graph(poly)
  T = Vector{Int32}[]
  for v in 1:length(graph)
    vnext = v == length(graph) ? 1 : v+1
    @assert vnext ∈ graph[v]
    push!(T,[v,vnext])
  end
  T
end

function Base.copy(poly::Polyhedron)
  vertices = get_vertex_coordinates(poly)
  graph = get_graph(poly)
  open = isopen(poly)
  data = get_metadata(poly)
  data = !isnothing(data) ? copy(data) : nothing
  Polyhedron(vertices,graph,open,data)
end

function simplexify_interior(p::Polygon)
  @assert !isopen(p)
  e_to_v = generate_facet_to_vertices(p)
  T = Vector{Int32}[]
  if length(e_to_v) > 0
    v0 = e_to_v[1][1]
    for verts in e_to_v
      if v0 ∉ verts
        push!(T,[v0,verts[1],verts[2]])
      end
    end
  end
  get_vertex_coordinates(p),T
end

"""
    simplexify_interior(p::Polyhedron)

  `simplex_interior` computes a simplex partition of the volume inside
  the Polyhedron `p`.
  It returns a vector of coordinates and an array of connectivitties.
"""
function simplexify_interior(poly::Polyhedron)
  !isopen(poly) || return simplexify_surface(poly)
  vstart = fill(UNSET,num_vertices(poly))
  stack = Int32[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    vstart[v] == UNSET || continue
    vstart[v] = v
    empty!(stack)
    push!(stack,v)
    while !isempty(stack)
      vcurrent = pop!(stack)
      for vneig in get_graph(poly)[vcurrent]
        if vstart[vneig] == UNSET
          vstart[vneig] = v
          push!(stack,vneig)
        end
      end
    end
  end
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  vertex_coordinates = get_vertex_coordinates(poly)
  T = Vector{Int32}[]
  X = vertex_coordinates
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        !isnothing(inext) || break
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        @assert vcurrent ≠ vnext
        vcurrent ≠ vnext || break
        if v ∉ (vstart[v],vcurrent,vnext)
          k = [vstart[v],v,vcurrent,vnext]
          push!(T,k)
        end
      end
    end
  end
  X,T
end

"""
    simplexify_surface(p::Polyhedron)

  `simplex_surface` computes a simplex partition of the surface bounding
  the Polyhedron `p`.
  It returns a vector of coordinates and an array of connectivitties.
"""
function simplexify_surface(poly::Polyhedron)
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  T = Vector{Int32}[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      vnext > 0 || continue
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        vnext > 0 || break
        if v ∉ (vcurrent,vnext)
          k = [v,vcurrent,vnext]
          push!(T,k)
        end
      end
    end
  end
  get_vertex_coordinates(poly),T
end
