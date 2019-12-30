module GenericBoundaryTriangulationsTests

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using FillArrays

import Gridap.Geometry: get_cell_coordinates
import Gridap.Geometry: get_cell_type
import Gridap.Geometry: get_reffes
import Gridap.Geometry: restrict
import Gridap.Geometry: get_normal_vector
import Gridap.Arrays: reindex

import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Fields: apply_lincomb
import Gridap.Fields: evaluate_field_array

abstract type BoundaryTriangulation{Dc,Dp} <: Triangulation{Dc,Dp} end

"""
"""
function get_volume_triangulation(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_face_to_cell(trian::BoundaryTriangulation)
  @abstractmethod
end

# and co

"""
"""
function get_face_to_cell_map(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function reindex(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function restrict(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function get_normal_vector(trian::BoundaryTriangulation)
  @abstractmethod
end

"""
"""
function BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
  GenericBoundaryTriangulation(model,face_to_mask)
end

function BoundaryTriangulation(model::DiscreteModel)
  topo = get_grid_topology(model)
  face_to_mask = collect(Bool,get_isboundary_face(topo))
  GenericBoundaryTriangulation(model,face_to_mask)
end


struct FaceToCellGlue{O} <: GridapType
  orientation::Val{O}
  face_to_cell::Vector{Int}
  face_to_lface::Vector{Int8}
  cell_to_lface_to_pindex::Table{Int8,Int32}
  ctype_to_lface_to_ftype::Vector{Vector{Int8}}
end

"""
Private fields
"""
struct GenericBoundaryTriangulation{Dc,Dp,Gf,Gc,O} <: BoundaryTriangulation{Dc,Dp}
  face_trian::Gf
  cell_trian::Gc
  glue::FaceToCellGlue{O}

  function GenericBoundaryTriangulation(
    face_trian::Triangulation,
    cell_trian::Triangulation,
    glue::FaceToCellGlue{O}) where O

    Dc = num_cell_dims(face_trian)
    Dp = num_point_dims(face_trian)
    Gf = typeof(face_trian)
    Gc = typeof(cell_trian)
    new{Dc,Dp,Gf,Gc,O}(face_trian, cell_trian, glue)
  end

end

function GenericBoundaryTriangulation(
  face_trian::Triangulation,
  cell_trian::Triangulation,
  topo::GridTopology,
  face_to_oldface::Vector{Int},
  icell_arround::Integer=1)


  D = num_cell_dims(model)
  oldface_to_cells = get_faces(topo,D-1,D)
  cell_to_oldfaces = get_faces(topo,D,D-1)
  cell_to_lface_to_pindex = Table(get_cell_permutations(topo,D-1))
  orientation = OrientationStyle(topo)

  oldface_to_cell = get_local_item(oldface_to_cells, icell_arround)
  oldface_to_lface = find_local_index(oldface_to_cell, cell_to_oldfaces)
  face_to_cell = collect(Int,reindex(oldface_to_cell, face_to_oldface))
  face_to_lface = collect(Int8,reindex(oldface_to_lface, face_to_oldface))

  f = (p)->fill(Int8(UNSET),num_faces(p,D-1))
  ctype_to_lface_to_ftype = map( f, get_reffes(cell_trian) )
  face_to_ftype = get_cell_type(face_trian)
  cell_to_ctype = get_cell_type(cell_trian)
  _fill_ctype_to_lface_to_ftype!(
    ctype_to_lface_to_ftype,
    face_to_cell,
    face_to_lface,
    face_to_ftype,
    cell_to_ctype)

  glue = FaceToCellGlue(
    orientation,
    face_to_cell,
    face_to_lface,
    cell_to_lface_to_pindex,
    ctype_to_lface_to_ftype)

  GenericBoundaryTriangulation(face_trian, cell_trian, glue)

end

function GenericBoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
  D = num_cell_dims(model)
  topo = get_grid_topology(model)
  oldface_grid = Grid(ReferenceFE{D-1},model)
  cell_grid = Grid(ReferenceFE{D},model)
  face_to_oldface = findall(get_isboundary_face(topo,1))
  face_grid = GridPortion(oldface_grid,face_to_oldface)
  GenericBoundaryTriangulation(face_grid,cell_grid,topo,face_to_oldface)
end

function _fill_ctype_to_lface_to_ftype!(
  ctype_to_lface_to_ftype,
  face_to_cell,
  face_to_lface,
  face_to_ftype,
  cell_to_ctype)
  for (face, cell) in enumerate(face_to_cell)
    ctype = cell_to_ctype[cell]
    lface = face_to_lface[face]
    ftype = face_to_ftype[face]
    ctype_to_lface_to_ftype[ctype][lface] = ftype
  end
end

# Triangulation Interface

function get_cell_coordinates(trian::GenericBoundaryTriangulation)
  get_cell_coordinates(trian.face_trian)
end

function get_reffes(trian::GenericBoundaryTriangulation)
    get_reffes(trian.face_trian)
end

function get_cell_type(trian::GenericBoundaryTriangulation)
  get_cell_type(trian.face_trian)
end

# Cell to face map

function get_face_to_cell_map(trian::GenericBoundaryTriangulation)
  face_to_fvertex_to_qcoors = FaceCellCoordinates(trian)
  f = (p)-> get_shapefuns(LagrangianRefFE(Float64,get_polytope(p),1))
  ftype_to_shapefuns = map( f, get_reffes(trian.face_trian) )
  face_to_ftype = collect(Int8, get_cell_type(trian.face_trian))
  face_to_shapefuns = CompressedArray(ftype_to_shapefuns, face_to_ftype)
  lincomb(face_to_shapefuns,face_to_fvertex_to_qcoors)
end

struct FaceCellCoordinates{D,T,O} <: AbstractVector{Vector{Point{D,T}}}
  glue::FaceToCellGlue{O}
  cell_to_ctype::Vector{Int8}
  ctype_to_lvertex_to_qcoords::Vector{Vector{Point{D,T}}}
  ctype_to_lface_to_lvertices::Vector{Vector{Vector{Int}}}
  ctype_to_lface_to_pindex_to_perm::Vector{Vector{Vector{Vector{Int}}}}
end

function FaceCellCoordinates(trian::GenericBoundaryTriangulation)
  d = num_cell_dims(trian)
  polytopes = map(get_polytope, get_reffes(trian.cell_trian))
  cell_to_ctype = collect(Int8,get_cell_type(trian.cell_trian))
  ctype_to_lvertex_to_qcoords = map(get_vertex_coordinates, polytopes)
  ctype_to_lface_to_lvertices = map((p)->get_faces(p,d,0), polytopes)
  ctype_to_lface_to_pindex_to_perm = map( (p)->get_face_vertex_permutations(p,d), polytopes)
  FaceCellCoordinates(
    trian.glue,
    cell_to_ctype,
    ctype_to_lvertex_to_qcoords,
    ctype_to_lface_to_lvertices,
    ctype_to_lface_to_pindex_to_perm)
end

function array_cache(a::FaceCellCoordinates{D,T,O}) where {D,T,O}
  if length(a.glue.face_to_cell) > 0
    cell = first(a.glue.face_to_cell)
    lface = first(a.glue.face_to_lface)
    ctype = a.cell_to_ctype[cell]
    lvertices = a.ctype_to_lface_to_lvertices[ctype][lface]
    c = zeros(Point{D,T},length(lvertices))
  else
    c = Point{D,T}[]
  end
  CachedArray(c)
end

function getindex!(cache,a::FaceCellCoordinates,face::Integer)
  cell = a.glue.face_to_cell[face]
  lface = a.glue.face_to_lface[face]
  ctype = a.cell_to_ctype[cell]
  cfvertex_to_lvertex = a.ctype_to_lface_to_lvertices[ctype][lface]
  p = a.glue.cell_to_lface_to_pindex.ptrs[cell]-1
  pindex = a.glue.cell_to_lface_to_pindex.data[p+lface]
  cfvertex_to_ffvertex = a.ctype_to_lface_to_pindex_to_perm[ctype][lface][pindex]
  lvertex_to_qcoods = a.ctype_to_lvertex_to_qcoords[ctype]
  nfvertices = length(cfvertex_to_lvertex)
  setsize!(cache,(nfvertices,))
  ffvertex_to_qcoords = cache.array
  for (cfvertex, ffvertex) in enumerate(cfvertex_to_ffvertex)
    lvertex = cfvertex_to_lvertex[cfvertex]
    qcoords = lvertex_to_qcoods[lvertex]
    ffvertex_to_qcoords[ffvertex] = qcoords
  end
  ffvertex_to_qcoords
end

function Base.getindex(a::FaceCellCoordinates,face::Integer)
  cache = array_cache(a)
  getindex!(cache,a,face)
end

Base.IndexStyle(::Type{FaceCellCoordinates{D,T,O}}) where {D,T,O} = IndexLinear()

Base.size(a::FaceCellCoordinates) = (length(a.glue.face_to_cell),)

# Value of face to cell map
# Not done for the moment, but the values can be precomputed, and only retrieved in the iteration

struct FaceToCellMapValue{T,V} <: AbstractVector{T}
  values::V
  function FaceToCellMapValue(ax::CompressedArray,b::FaceCellCoordinates)
    values = apply_lincomb_default(ax,b)
    V = typeof(values)
    T = eltype(V)
    new{T,V}(values)
  end
end

Base.size(v::FaceToCellMapValue) = size(v.values)

Base.IndexStyle(::Type{FaceToCellMapValue{T,V}}) where {T,V} = IndexStyle(V)

Base.getindex(v::FaceToCellMapValue,i::Integer) = v.values[i]

array_cache(v::FaceToCellMapValue) = array_cache(v.values)

getindex!(cache,v::FaceToCellMapValue,i::Integer) = getindex!(cache,v.values,i)

function apply_lincomb(ax::CompressedArray,b::FaceCellCoordinates)
  FaceToCellMapValue(ax,b)
end

struct CellShapeFunsAtFaces{T} <: AbstractVector{T}
end

function evaluate_field_array(g::CompressedArray{<:Field},fx::FaceToCellMapValue)
  apply(g,fx)
end

# Tests

domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
model = UnstructuredDiscreteModel(model)
topo = get_grid_topology(model)

btrian = BoundaryTriangulation(model)
test_triangulation(btrian)

btrian = GenericBoundaryTriangulation(model,get_isboundary_face(topo,1))
test_triangulation(btrian)

face_to_fvertex_to_qcoors = FaceCellCoordinates(btrian)
r =  Vector{Point{2,Float64}}[
  [(0.0,0.0),(1.0,0.0)],[(0.0,0.0),(0.0,1.0)],
  [(0.0,0.0),(1.0,0.0)],[(1.0,0.0),(1.0,1.0)],
  [(0.0,1.0),(1.0,1.0)],[(0.0,0.0),(0.0,1.0)],
  [(0.0,1.0),(1.0,1.0)],[(1.0,0.0),(1.0,1.0)]]
test_array(face_to_fvertex_to_qcoors,r)

s2q = get_face_to_cell_map(btrian)

s = CompressedArray([Point{1,Float64}[(0.25,),(0.75,)]],get_cell_type(btrian))

q = evaluate(s2q,s)

q2x = reindex(get_cell_map(get_grid(model)), btrian.glue.face_to_cell)

s2x = compose(q2x,s2q)

x = evaluate(s2x,s)

display(s)
display(q)
display(x)
display((eltype(get_cell_map(get_grid(model)))))


end # module
