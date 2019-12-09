module DiscreteModelsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock

model = DiscreteModelMock()
test_discrete_model(model)
@test is_oriented(model) == false

@test num_dims(model) == 2
@test num_cell_dims(model) == 2
@test num_point_dims(model) == 2
@test num_faces(model) == 9 + 13 + 5
@test num_faces(model,0) == 9
@test num_faces(model,1) == 13
@test num_faces(model,2) == 5
@test num_cells(model) == 5
@test num_vertices(model) == 9
@test num_edges(model) == 13
@test num_facets(model) == 13
@test num_nodes(model) == 9
@test get_vertex_coordinates(model) == get_node_coordinates(model)
@test get_isboundary_face(model) == Bool[1,1,1,1,0,1,1,1,1,1,0,1,0,1,0,1,0,1,1,0,1,1,1,1,1,1,1]
@test get_isboundary_node(model) == Bool[1,1,1,1,0,1,1,1,1]
@test get_dimranges(model) == [1:9, 10:22, 23:27]
@test get_offsets(model) == [0, 9, 22]
@test get_offset(model,0) == 0
@test get_offset(model,1) == 9
@test get_offset(model,2) == 22
@test get_facedims(model) == Int8[0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2]
@test get_cell_faces(model) == [
  [1,2,4,5,10,11,12,13,23], [2,3,5,14,13,15,24], [3,6,5,16,15,17,25],
  [4,5,7,8,11,18,19,20,26], [5,6,8,9,17,21,20,22,27]]

grid = UnstructuredGrid(ReferenceFE{0},model)
@test num_cells(grid) == num_vertices(model)

grid = ConformingTriangulation(ReferenceFE{1},model)
@test num_cells(grid) == num_edges(model)

grid = ConformingTriangulation(ReferenceFE{2},model)
@test num_cells(grid) == num_cells(model)


using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry: _get_cell_data
import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!


function replace_reffes(model::DiscreteModel,reffes::Vector{<:NodalReferenceFE})

  face_to_own_nodes, n_nodes = _generate_face_to_own_nodes(
    model,reffes,ConformityStyle(model))

  cell_to_nodes = _generate_cell_to_nodes(face_to_own_nodes,model,reffes,)

end

function _generate_face_to_own_nodes(model,reffes,conf)
  @notimplemented
end

function _generate_face_to_own_nodes(model,reffes,conf::RegularConformity)

  D = num_cell_dims(model)
  cell_to_ctype = get_cell_type(model)

  d_to_ctype_to_ldface_to_own_lnodes = [
    [ get_face_own_nodes(reffe,d) for reffe in reffes] for d in 0:D]

  d_to_dface_to_cells = [ get_faces(model,d,D) for d in 0:D]
  d_to_cell_to_dfaces = [ get_faces(model,D,d) for d in 0:D]

  d_to_offset = get_offsets(model)

  n_faces = num_faces(model)

  face_to_own_nodes, n_nodes =  _generate_face_to_own_dofs(
    n_faces,
    cell_to_ctype,
    d_to_cell_to_dfaces,
    d_to_dface_to_cells,
    d_to_offset,
    d_to_ctype_to_ldface_to_own_lnodes)

  (face_to_own_nodes, n_nodes)

end

function _generate_face_to_own_dofs(
  n_faces,
  cell_to_ctype,
  d_to_cell_to_dfaces::Vector{Table{T,P}},
  d_to_dface_to_cells::Vector{Table{T,P}},
  d_to_offset,
  d_to_ctype_to_ldface_to_own_ldofs) where {T,P}

  face_to_own_dofs_ptrs = zeros(P,n_faces+1)

  D = length(d_to_offset)-1

  for d in 0:D
    cell_to_dfaces = d_to_cell_to_dfaces[d+1]
    dface_to_cells = d_to_dface_to_cells[d+1]
    offset = d_to_offset[d+1]
    ctype_to_ldface_to_own_ldofs = d_to_ctype_to_ldface_to_own_ldofs[d+1]
    ctype_to_ldface_to_num_own_ldofs = map( (x) -> length.(x) ,ctype_to_ldface_to_own_ldofs)
    icell = 1
    dface_to_cell_owner = get_local_item(dface_to_cells,icell)
    dface_to_ldface = find_local_index(dface_to_cell_owner,cell_to_dfaces)

    if any( ctype_to_ldface_to_num_own_ldofs .!= 0)
      _generate_face_to_own_dofs_count_d!(
        face_to_own_dofs_ptrs,
        offset,
        cell_to_ctype,
        dface_to_cell_owner,
        dface_to_ldface,
        ctype_to_ldface_to_num_own_ldofs)
    end
  end

  length_to_ptrs!(face_to_own_dofs_ptrs)
  n_dofs = face_to_own_dofs_ptrs[end]-1
  face_to_own_dofs_data = collect(T(1):T(n_dofs))

  face_to_own_dofs = Table(face_to_own_dofs_data,face_to_own_dofs_ptrs)
  (face_to_own_dofs, n_dofs)
end

function  _generate_face_to_own_dofs_count_d!(
  face_to_own_dofs_ptrs,
  offset,
  cell_to_ctype,
  dface_to_cell_owner,
  dface_to_ldface,
  ctype_to_ldface_to_num_own_ldofs)

  n_dfaces = length(dface_to_ldface)
  for dface in 1:n_dfaces
    cell = dface_to_cell_owner[dface]
    ldface = dface_to_ldface[dface]
    ctype = cell_to_ctype[cell]
    n_own_ldofs = ctype_to_ldface_to_num_own_ldofs[ctype][ldface]
    face = dface + offset
    face_to_own_dofs_ptrs[face+1] = n_own_ldofs
  end
end

function _generate_cell_to_nodes(
  face_to_own_nodes, model, reffes, orientation, conformity)
  @notimplemented
end

function _generate_cell_to_nodes(
  face_to_own_nodes, model, reffes, orientation::Val{true}, ::RegularConformity)

  cell_to_faces = get_cell_faces(model)
  cell_to_ctype = get_cell_type(model)
  ctype_lface_to_own_lnodes = [get_face_own_nodes(reffe) for reffe in reffes]
  ctype_to_num_dofs = map(num_dofs,reffes)

  cell_to_nodes = CellDofsOriented(
    cell_to_faces,
    cell_to_ctype,
    ctype_lface_to_own_lnodes,
    ctype_to_num_dofs,
    face_to_own_nodes)

  Table(cell_to_nodes)

end

struct CellDofsOriented{T,P,V<:AbstractVector} <: AbstractVector{Vector{T}}
  cell_to_faces::Table{T,P}
  cell_to_ctype::V
  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_to_num_dofs::Vector{Int}
  face_to_own_dofs::Table{T,P}
end

Base.size(a::CellDofsOriented) = size(a.cell_to_faces)

Base.IndexStyle(::Type{<:CellDofsOriented}) = IndexStyle(Table)

function array_cache(a::CellDofsOriented)
  n_dofs = testitem(a.ctype_to_num_dofs)
  T = eltype(eltype(a))
  v = zeros(T,n_dofs)
  CachedArray(v)
end

function getindex!(cache,a::CellDofsOriented,cell::Integer)
  ctype = a.cell_to_ctype[cell]
  n_dofs = a.ctype_to_num_dofs[ctype]
  setsize!(cache,(n_dofs,))
  dofs = cache.array
  lface_to_own_ldofs = a.ctype_to_lface_to_own_ldofs[ctype]
  p = a.cell_to_faces.ptrs[cell]-1
  for (lface, own_ldofs) in enumerate(lface_to_own_ldofs)
    face = a.cell_to_faces.data[p+lface]
    q = a.face_to_own_dofs.ptrs[face]-1
    for (i,ldof) in enumerate(own_ldofs)
      dof = a.face_to_own_dofs.data[q+i]
      dofs[ldof] = dof
    end
  end
  dofs
end

function Base.getindex(a::CellDofsOriented,cell::Integer)
  cache = array_cache(a)
  getindex!(cache,a,cell)
end

function _generate_node_to_coords(n_nodes, cell_to_lnode_to_node, model,reffes,conformity)
  @notimplemented
end

function _generate_node_to_coords(n_nodes, cell_to_lnode_to_node, model, reffes, ::RegularConformity)

  cell_to_map = get_cell_map(model)
  ctype_to_lnode_to_lcoords = map(get_node_coordinates,reffes)
  cell_to_ctype = get_cell_type(model)
  cell_to_lnode_to_lcoords = _get_cell_data(ctype_to_lnode_to_lcoords, cell_to_ctype)
  cell_to_lnode_to_coords = evaluate(cell_to_map,cell_to_lnode_to_lcoords)

  T = eltype(eltype(get_node_coordinates(model)))
  D = num_point_dims(model)
  node_to_coords = zeros(Point{D,T},n_nodes)

  cache = array_cache(cell_to_lnode_to_coords)

  _fill_node_to_coords!(
    node_to_coords,
    cell_to_lnode_to_node,
    cell_to_lnode_to_coords,
    cache)

  node_to_coords

end

function  _fill_node_to_coords!(
  node_to_coords,
  cell_to_lnode_to_node::Table,
  cell_to_lnode_to_coords,
  cache)
  for cell in 1:length(cell_to_lnode_to_node)
    p = cell_to_lnode_to_node.ptrs[cell]-1
    lnode_to_coords = getindex!(cache,cell_to_lnode_to_coords,cell)
    for (lnode, coords) in enumerate(lnode_to_coords)
      node = cell_to_lnode_to_node.data[p+lnode]
      node_to_coords[node] = coords
    end
  end
end

domain = (0,1,0,1)
partition = (2,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)
order = 3

reffes = [ LagrangianRefFE(Float64,get_polytope(reffe),order) for reffe in get_reffes(model)]

face_to_own_nodes, n_nodes = _generate_face_to_own_nodes(model,reffes,ConformityStyle(model))

@show n_nodes
display(face_to_own_nodes)

cell_to_nodes = _generate_cell_to_nodes(
  face_to_own_nodes, model, reffes,
  OrientationStyle(model), ConformityStyle(model))

display(cell_to_nodes)

node_to_coords = _generate_node_to_coords(n_nodes, cell_to_nodes, model, reffes, ConformityStyle(model))

display(node_to_coords)

end # module
