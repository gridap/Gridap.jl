module UnstructuredGridTopologiesTests

using Test
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: GridTopologyMock

import Gridap.Geometry: num_faces
import Gridap.Geometry: get_faces
import Gridap.Geometry: get_cell_type
import Gridap.Geometry: get_polytopes
import Gridap.Geometry: get_vertex_coordinates

using Gridap.Geometry: generate_cells_around
using Gridap.Geometry: generate_cell_to_faces
using Gridap.Geometry: generate_face_to_vertices

struct UnstructuredGridTopology{Dc,Dp,T} <: GridTopology{Dc,Dp}
  vertex_coordinates::Vector{Point{Dp,T}}
  n_m_to_nface_to_mfaces::Matrix{Table{Int,Int32}}
  cell_type::Vector{Int8}
  polytopes::Vector{Polytope{Dc}}
end

# Constructors

function UnstructuredGridTopology(
  vertex_coordinates::Vector{<:Point},
  cell_vertices::Table,
  cell_type::Vector{<:Integer},
  polytopes::Vector{<:Polytope})

  D = num_dims(first(polytopes))
  n = D+1
  n_m_to_nface_to_mfaces = Matrix{Table{Int,Int32}}(undef,n,n)
  n_m_to_nface_to_mfaces[D+1,0+1] = cell_vertices
  vertex_cells = generate_cells_around(cell_vertices,length(vertex_coordinates))
  n_m_to_nface_to_mfaces[0+1,D+1] = vertex_cells

  P = eltype(vertex_coordinates)
  Dp = length(P)
  T = eltype(P)

  UnstructuredGridTopology{D,Dp,T}(
    vertex_coordinates,
    n_m_to_nface_to_mfaces,
    cell_type,
    polytopes)

end

# Needed, do not remove
function num_faces(g::UnstructuredGridTopology,d::Integer)
  if d == 0
    length(g.vertex_coordinates)
  elseif d == num_cell_dims(g)
    length(g.cell_type)
  else
    D = num_cell_dims(g)
    face_to_cells = get_faces(g,d,D)
    length(face_to_cells)
  end
end

# Implementation of abstract API

get_vertex_coordinates(g::UnstructuredGridTopology) = g.vertex_coordinates

get_cell_type(g::UnstructuredGridTopology) = g.cell_type

get_polytopes(g::UnstructuredGridTopology) = g.polytopes

function get_faces(g::UnstructuredGridTopology, dimfrom::Integer, dimto::Integer)
  _setup_faces!(g,dimfrom,dimto)
  g.n_m_to_nface_to_mfaces[dimfrom+1,dimto+1]
end

function _setup_faces!(g,dimfrom,dimto)
  if isassigned(g.n_m_to_nface_to_mfaces,dimfrom+1,dimto+1)
    return nothing
  end
  D = num_cell_dims(g)
  if dimfrom==0
    if dimto==0
      _setup_face_to_face!(g,dimfrom)
    elseif dimto==D
      nothing
    else
      _setup_face_to_vertices!(g,dimto)
    end
  elseif dimfrom == D
    if dimto==0
      nothing
    elseif dimto==D
      _setup_face_to_face!(g,dimfrom)
    else
      _setup_cell_to_faces!(g,dimto)
    end
  else
    if dimto==0
      _setup_face_to_vertices!(g,dimfrom)
    elseif dimto==D
      _setup_cell_to_faces!(g,dimfrom)
    elseif dimto==dimfrom
      _setup_face_to_face!(g,dimfrom)
    elseif dimfrom > dimto
      _setup_nface_to_mface!(g,dimfrom,dimto)
    else
      _setup_nface_to_mface!(g,dimto,dimfrom)
    end
  end
  nothing
end

function _setup_cell_to_faces!(model,dimto)

  D = num_cell_dims(model)

  if isassigned(model.n_m_to_nface_to_mfaces,D+1,dimto+1)
    return
  end
  if isassigned(model.n_m_to_nface_to_mfaces,dimto+1,0+1)
    @notimplemented

  else

    cell_to_vertices = get_faces(model,D,0)
    vertex_to_cells = get_faces(model,0,D)
    cell_to_cell_type = get_cell_type(model)
    polytopes = get_polytopes(model)
    cell_type_to_lface_to_lvertices = map( (p)->get_faces(p,dimto,0), polytopes )

    cell_to_faces = generate_cell_to_faces(
        cell_to_vertices,
        cell_type_to_lface_to_lvertices,
        cell_to_cell_type,
        vertex_to_cells)

    faces_to_cells = generate_cells_around(cell_to_faces)

    model.n_m_to_nface_to_mfaces[D+1,dimto+1] = cell_to_faces
    model.n_m_to_nface_to_mfaces[dimto+1,D+1] = faces_to_cells

  end

  nothing

end

function _setup_face_to_vertices!(model,dimfrom)

  if isassigned(model.n_m_to_nface_to_mfaces,dimfrom+1,0+1)
    return
  end

  D = num_cell_dims(model)

  cell_to_vertices = get_faces(model,D,0)
  cell_to_faces = get_faces(model,D,dimfrom)
  cell_to_ctype = get_cell_type(model)
  polytopes = get_polytopes(model)
  ctype_to_lface_to_lvertices = map( (p)->get_faces(p,dimfrom,0), polytopes )

  nfaces = num_faces(model,dimfrom)

  face_to_vertices = generate_face_to_vertices(
    cell_to_vertices,
    cell_to_faces,
    cell_to_ctype,
    ctype_to_lface_to_lvertices,
    nfaces)

  nvertices = num_faces(model,0)

  vertex_to_faces = generate_cells_around(face_to_vertices,nvertices)

  model.n_m_to_nface_to_mfaces[dimfrom+1,0+1] = face_to_vertices
  model.n_m_to_nface_to_mfaces[0+1,dimfrom+1] = vertex_to_faces

  nothing

end

function _setup_face_to_face!(model,d)

  if isassigned(model.n_m_to_nface_to_mfaces,d+1,d+1)
    return
  end

  nfaces = num_faces(model,d)
  id = identity_table(Int,Int32,nfaces)
  model.n_m_to_nface_to_mfaces[d+1,d+1] = id

  return

end

function _setup_nface_to_mface!(model,n,m)
  @assert n > m

  if isassigned(model.n_m_to_nface_to_mfaces,n+1,m+1)
    return
  end

  nface_to_vertices = get_faces(model,n,0)
  vertex_to_mfaces = get_faces(model,0,m)
  nface_to_nftype = get_face_type(model,n)

  polytopes = get_reffaces(Polytope{n},model)
  nftype_to_lmface_to_lvertices = map( (p)->get_faces(p,m,0), polytopes)

  nface_to_mfaces = find_cell_to_faces(
    nface_to_vertices,
    nftype_to_lmface_to_lvertices,
    nface_to_nftype,
    vertex_to_mfaces)

  num_mfaces = num_faces(model,m)
  mface_to_nfaces = generate_cells_around(nface_to_mfaces,num_mfaces)

  model.n_m_to_nface_to_mfaces[n+1,m+1] = nface_to_mfaces
  model.n_m_to_nface_to_mfaces[m+1,n+1] = mface_to_nfaces

  nothing

end

m = GridTopologyMock()

g = UnstructuredGridTopology(
  get_vertex_coordinates(m),
  get_cell_vertices(m),
  get_cell_type(m),
  get_polytopes(m))

test_grid_topology(g)

@test num_faces(g,0) == num_faces(m,0)
@test num_faces(g,1) == num_faces(m,1)
@test num_faces(g,2) == num_faces(m,2)
@test get_reffaces(Polytope{1},g) == get_reffaces(Polytope{1},m)
@test get_reffaces(Polytope{0},g) == get_reffaces(Polytope{0},m)
@test get_reffaces(Polytope{2},g) == get_reffaces(Polytope{2},m)
@test get_face_type(g,1) == get_face_type(m,1)
@test get_isboundary_face(g,0) == get_isboundary_face(m,0)
@test get_isboundary_face(g,1) == get_isboundary_face(m,1)
@test get_isboundary_face(g,2) == get_isboundary_face(m,2)
@test get_faces(g,2,0) == get_faces(m,2,0)
@test get_faces(g,2,1) == get_faces(m,2,1)
@test get_faces(g,1,2) == get_faces(m,1,2)
@test get_faces(g,2,1) == get_faces(m,2,1)
@test get_faces(g,0,1) == get_faces(m,0,1)
@test get_faces(g,1,0) == get_faces(m,1,0)
@test get_faces(g,1,1) == get_faces(m,1,1)
@test get_vertex_coordinates(g) == get_vertex_coordinates(m)           

end # module
