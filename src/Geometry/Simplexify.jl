module Simplexify

using Gridap
using Gridap.Helpers
using UnstructuredGrids.Kernels: refine_grid_connectivity
using UnstructuredGrids.Kernels: generate_tface_to_face
using Gridap.DiscreteModels: DiscreteModelFromData

export simplexify

function simplexify(grid::Grid{D,Z}) where {D,Z}
  ugrid = UnstructuredGrid(grid)

  ltcell_to_lpoints = _generate_ltcell_to_lpoints(celltypes(ugrid))

  tcells_data, tcells_ptrs = refine_grid_connectivity(
    ugrid.cells_data, ugrid.cells_ptrs, ltcell_to_lpoints)

  order = 1
  _check_order(order,cellorders(ugrid))

  ntcells = length(tcells_ptrs) - 1
  ex = tuple(fill(TET_AXIS,Z)...)
  _ct = ConstantCellValue(ex,ntcells)
  _co = ConstantCellValue(order,ntcells)

  UnstructuredGrid(
    points(ugrid),
    tcells_data,
    tcells_ptrs,
    _ct,
    _co)

end

function simplexify(model::DiscreteModel{D}) where D

  grid = Grid(model)

  graph = FullGridGraph(model)

  facelabels = FaceLabels(model)

  tgrid = simplexify(grid)

  tgraph = FullGridGraph(tgrid)

  tgrids = _generate_tgrids(tgrid,tgraph)

  ltcell_to_lnodes = _generate_ltcell_to_lpoints(celltypes(grid))

  tfacelabels = _generate_tfacelabels(
    tgrids,grid,facelabels,tgraph,graph,ltcell_to_lnodes)

  DiscreteModelFromData(tgrids,tgraph,tfacelabels)

end

# Helpers

const UNSET_ID = 0

function _generate_tfacelabels(
  tgrids,grid,facelabels,tgraph,graph,ltcell_to_lnodes)

  dim_to_tface_to_label = [ fill(UNSET_ID,ncells(g)) for g in tgrids ]

  dim_to_face_to_label = facelabels.dim_to_nface_to_label

  tgrid = tgrids[end]
  tpolytope = CellPolytopes(tgrid).value
  polytope = CellPolytopes(grid).value

  _fill_dim_to_tface_to_label!(
    dim_to_tface_to_label,
    dim_to_face_to_label,
    tpolytope,
    polytope,
    tgraph,
    graph,
    ltcell_to_lnodes)

  FaceLabels(
    dim_to_tface_to_label, facelabels.tag_to_labels, facelabels.tag_to_name)

end

function _fill_dim_to_tface_to_label!(
  dim_to_tface_to_label,
  dim_to_face_to_label,
  tpolytope,
  polytope,
  tgraph,
  graph,
  ltcell_to_lnodes)

  D = length(dim_to_face_to_label)-1
  d = 0
  dim_to_tface_to_label[d+1] = dim_to_face_to_label[d+1]

  trefcell = RefCell(tpolytope)
  refcell = RefCell(polytope)

  for d in 1:(D-1)

    cell_to_faces = connections(graph,D,d)
    cell_to_faces_data, cell_to_faces_ptrs = compress(cell_to_faces)

    tcell_to_tfaces = connections(tgraph,D,d)

    ntfaces = maximum(tcell_to_tfaces.data)

    ltface_to_ltnodes = trefcell.faces[d+1]
    lface_to_lnodes = refcell.faces[d+1]

    tface_to_face = generate_tface_to_face(
      cell_to_faces_data,
      cell_to_faces_ptrs,
      tcell_to_tfaces.data,
      tcell_to_tfaces.ptrs,
      ltcell_to_lnodes,
      ltface_to_ltnodes,
      lface_to_lnodes,
      ntfaces)

    _update_labels!(
      dim_to_tface_to_label[d+1],dim_to_face_to_label[d+1],tface_to_face)

  end

  d = D
  ncells = length(dim_to_face_to_label[d+1])
  nltcells = length(ltcell_to_lnodes)
  tcell_to_cell = _generate_tcell_to_cell(ncells,nltcells)
  _update_labels!(
    dim_to_tface_to_label[d+1],dim_to_face_to_label[d+1],tcell_to_cell)

  for d = 1:(D-1)

    for j in (d+1):D

      dface_to_jfaces = connections(tgraph,d,j)
      dface_to_label = dim_to_tface_to_label[d+1]
      jface_to_label = dim_to_tface_to_label[j+1]
      _fix_dface_to_label!(dface_to_label,jface_to_label,dface_to_jfaces)

    end

  end

end

function _fix_dface_to_label!(dface_to_label,jface_to_label,dface_to_jfaces)

  ndfaces = length(dface_to_label)
  @assert ndfaces == length(dface_to_jfaces)

  for dface in 1:ndfaces

    dlabel = dface_to_label[dface]
    if dlabel != UNSET_ID
      continue
    end

    jfaces = dface_to_jfaces[dface]
    for jface in jfaces
      jlabel = jface_to_label[jface]
      if jlabel != UNSET_ID
        dface_to_label[dface] = jlabel
        break
      end
    end

  end

end

function _update_labels!(tface_to_label,face_to_label,tface_to_face)
  for tface in 1:length(tface_to_label)
    face = tface_to_face[tface]
    if face != 0
      tface_to_label[tface] = face_to_label[face]
    end
  end
end

function _generate_tcell_to_cell(ncells,nltcells)
  ntcells = ncells * nltcells
  tcell_to_cell = Vector{Int}(undef,ntcells)
  tcell = 1
  for cell in 1:ncells
    for ltcell in 1:nltcells
      tcell_to_cell[tcell] = cell
      tcell += 1
    end
  end
  tcell_to_cell
end

function _generate_tgrids(tgrid::Grid{D},tgraph) where D
  tgrids = Grid[]
  for d in 0:(D-1)
    tface_to_vertices = connections(tgraph,d,0)
    tgrid_d = _grid_of_dim(points(tgrid),tface_to_vertices,d)
    push!(tgrids,tgrid_d)
  end
  push!(tgrids,tgrid)
  tgrids
end

function _grid_of_dim(coords,face_to_vertices,Z)

  nfaces = length(face_to_vertices)
  fcode = tuple([TET_AXIS for i in 1:Z]...)

  order = 1

  _points = coords
  _cells_data = face_to_vertices.data
  _cells_ptrs = face_to_vertices.ptrs
  _ctypes = ConstantCellValue(fcode,nfaces)
  _corders = ConstantCellValue(order,nfaces)

  UnstructuredGrid(
    _points, _cells_data, _cells_ptrs, _ctypes, _corders)

end

function _check_order(order,co)
  @notimplemented
end

function _check_order(order,co::ConstantCellValue)
  @notimplementedif co.value != order
end

function _generate_ltcell_to_lpoints(ct)
  @notimplemented
end

function _generate_ltcell_to_lpoints(ct::ConstantCellValue)
  extrusion = ct.value
  if extrusion == (HEX_AXIS, HEX_AXIS)
    ltcell_to_lpoints = [[1,2,3],[4,3,2]]
    return ltcell_to_lpoints
  elseif extrusion == (HEX_AXIS, HEX_AXIS, HEX_AXIS)
    ltcell_to_lpoints = [
      [7,3,2,1], [7,5,2,1], [7,4,3,2], [7,4,8,2], [7,6,5,2], [7,6,8,2]]
  else
    @notimplemented
  end
end

end # module
