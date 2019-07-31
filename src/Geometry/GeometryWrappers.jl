module GeometryWrappers #@fverdugo rename as RefCells ??

# Dependencies of this module

using Gridap
using Gridap.Helpers
using UnstructuredGrids.Core: VERTEX
using UnstructuredGrids.Kernels: UNSET
using UnstructuredGrids.Kernels: generate_data_and_ptrs

# Functionality provided by this module

export RefCell
import UnstructuredGrids.Core: RefCell
import Gridap: Grid

"""
Construct a RefCell from a Polytope
"""
RefCell(polytope::Polytope) = _ref_cell_from_polytope(polytope)

"""
Construct a Grid from a polytope
"""
function Grid(polytope::Polytope{D},dim::Int) where D
  if dim < D
    _grid_d(polytope,dim)
  elseif dim==D
    _grid_D(polytope)
  else
    @unreachable
  end
end

# Helpers

function _grid_D(polytope::Polytope{D}) where D

  order = 1
  reffe = LagrangianRefFE{D,Float64}(polytope,order)
  points = reffe.dofbasis.nodes

  cell_to_nodes = [ [i for i in 1:length(points)], ]

  cells_data,cells_ptrs = generate_data_and_ptrs(cell_to_nodes)
  extrusion = polytope.extrusion.array.data

  l = 1
  co = ConstantCellValue(order,l)
  ct = ConstantCellValue(extrusion,l)

  UnstructuredGrid(points,cells_data,cells_ptrs,ct,co)

end

function _grid_d(polytope::Polytope{D},dim::Int) where D
  @assert dim < D

  order = 1
  reffe = LagrangianRefFE{D,Float64}(polytope,order)
  points = reffe.dofbasis.nodes

  dim_to_jface_to_vertices, dim_to_jface_to_code = _faces(polytope)
  jface_to_vertices = dim_to_jface_to_vertices[dim+1]
  jface_to_code = dim_to_jface_to_code[dim+1]

  njfaces = length(jface_to_code)
  @assert njfaces > 0
  code1 = jface_to_code[1]
  # TODO for the moment, we assume all faces of the polytope have the same
  # extrusion code
  @notimplementedif any([ code1 != code for code in jface_to_code ])

  cells_data, cells_ptrs = generate_data_and_ptrs(jface_to_vertices)
  code = (code1...,)
  ctypes = ConstantCellValue(code,njfaces)
  order = 1
  orders = ConstantCellValue(order,njfaces)

  UnstructuredGrid(points,cells_data,cells_ptrs,ctypes,orders)

end



function _ref_cell_from_polytope(polytope::Polytope{D}) where D

  dim_to_jface_to_vertices, dim_to_jface_to_code = _faces(polytope)

  vertex_to_coords = _coordinates(polytope)

  dim_to_jface_to_ftype, dim_to_ftype_to_code = _ftypes(dim_to_jface_to_code)

  dim_to_ftype_to_refface = _reffaces(dim_to_ftype_to_code)

  _vtkid, _vtknodes = _vtkinfo(polytope)

  RefCell(
    ndims = D,
    faces = dim_to_jface_to_vertices,
    facetypes = dim_to_jface_to_ftype,
    reffaces = dim_to_ftype_to_refface,
    coordinates = vertex_to_coords,
    vtkid = _vtkid,
    vtknodes = _vtknodes)

end

function _faces(polytope)

  nface_to_dim, nface_to_vertices, nface_to_code = (
    _extract_nface_to_info(polytope) )

 dim_to_jface_to_vertices, dim_to_jface_to_code = (
    _sort_nfaces_to_vertices_per_dim(
      nface_to_dim, nface_to_vertices, nface_to_code) )

 (dim_to_jface_to_vertices, dim_to_jface_to_code)
end

function _reffaces(dim_to_ftype_to_code)

  [ [ _code_to_refface(code) for code in ftype_to_code ]
   for ftype_to_code in dim_to_ftype_to_code ]

end

function _code_to_refface(code)
  if length(code) == 0
    return VERTEX
  else
    t = tuple(code...)
    p = Polytope(t)
    return RefCell(p)
  end
end

function _ftypes(dim_to_jface_to_code)
  dim_to_jface_to_ftype = Vector{Int}[]
  dim_to_ftype_to_code = Vector{Vector{Int}}[]
  for d in 1:length(dim_to_jface_to_code)
    jface_to_code = dim_to_jface_to_code[d]
    jface_to_ftype, ftype_to_code = _find_unique_codes(jface_to_code)
    push!(dim_to_jface_to_ftype,jface_to_ftype)
    push!(dim_to_ftype_to_code,ftype_to_code)
  end
  (dim_to_jface_to_ftype, dim_to_ftype_to_code)
end

function _coordinates(polytope)
  num_nfaces = length(polytope.nf_dim)
  nvertices = 0
  for nface in 1:num_nfaces
    d = length(polytope.nf_dim[nface])-1
    if d == 0
      nvertices += 1
    end
  end
  dim = length(polytope.extrusion)
  coords = Array{Float64,2}(undef,(dim,nvertices))
  x = Vector{Float64}(undef,dim)
  for nface in 1:num_nfaces
    d = length(polytope.nf_dim[nface])-1
    if d == 0
      x .= polytope.nfaces[nface].anchor.array
      coords[:,nface] .= x
    end
  end
  coords
end


# TODO use this once NodesArray works for general polytopes
#function _coordinates(polytope::Polytope{D}) where D
#  orders = fill(1,D)
#  na = NodesArray(polytope,orders)
#  x = na.coordinates
#  nnodes = length(x)
#  coords = Array{Float64,2}(undef,(D,nnodes))
#  for i in 1:nnodes
#    for d in 1:D
#      coords[d,i] = x[i][d]
#    end
#  end
#  coords
#end

function _vtknodes(code)
  h = HEX_AXIS
  t = TET_AXIS
  if code == [] ; [1,]
  elseif code[2:end] == [] ; [1,2]
  elseif code[2:end] == [t] ; [1,2,3]
  elseif code[2:end] == [h] ; [1,2,4,3]
  elseif code[2:end] == [t,t] ; [1,2,3,4]
  elseif code[2:end] == [h,h] ; [1,2,4,3,5,6,8,7]
  else; Int[]
  end
end

function _vtkid(code)
  h = HEX_AXIS
  t = TET_AXIS
  if code == []; 1
  elseif code[2:end] == []; 3
  elseif code[2:end] == [t]; 5
  elseif code[2:end] == [h]; 9
  elseif code[2:end] == [t,t]; 10
  elseif code[2:end] == [h,h]; 12
  else; UNSET
  end
end

function _vtkinfo(polytope)
  code = [i for i in polytope.extrusion.array]
  (_vtkid(code), _vtknodes(code))
end

function _extract_nface_to_info(polytope)
  num_nfaces = length(polytope.nf_dim)
  nface_to_dim = zeros(Int,num_nfaces)
  nface_to_vertices = Vector{Vector{Int}}(undef,num_nfaces)
  nface_to_code = Vector{Vector{Int}}(undef,num_nfaces)
  vdim = 0
  for nface in 1:num_nfaces
    nface_to_dim[nface] = length(polytope.nf_dim[nface])-1
    v = polytope.nf_dim[nface][vdim+1]
    nface_to_vertices[nface] = polytope.nf_nfs[nface][v]
    nface_to_code[nface] = _nface_code(polytope.nfaces[nface])
  end
  (nface_to_dim, nface_to_vertices, nface_to_code)
end

function _nface_code(nface)
  [ i for i in nface.extrusion if i != 0 ]
end

function _sort_nfaces_to_vertices_per_dim(
  nface_to_dim, nface_to_vertices, nface_to_code)
  num_nfaces = length(nface_to_dim)
  ndims = maximum(nface_to_dim)
  dim_to_jface_to_vertices = Vector{Vector{Vector{Int}}}(undef,ndims)
  dim_to_jface_to_code = Vector{Vector{Vector{Int}}}(undef,ndims)
  vdim = 0
  for d in 0:ndims-1
    jface_to_vertices = Vector{Vector{Int}}(undef,0)
    jface_to_code = Vector{Vector{Int}}(undef,0)
    for nface in 1:num_nfaces
      if nface_to_dim[nface] == d
        push!(jface_to_vertices,nface_to_vertices[nface])
        push!(jface_to_code, nface_to_code[nface])
      end
    end
    dim_to_jface_to_vertices[d+1] = jface_to_vertices
    dim_to_jface_to_code[d+1] = jface_to_code
  end
  dim_to_jface_to_vertices, dim_to_jface_to_code
end

function _find_unique_codes(jface_to_code)
  ftype_to_code = unique(jface_to_code)
  _jface_to_ftype = indexin(jface_to_code,ftype_to_code)
  jface_to_ftype = _remove_nothing(_jface_to_ftype)
  (jface_to_ftype, ftype_to_code)
end

function _remove_nothing(a)
  b = Vector{Int}(undef,size(a))
  b .= a
end

end # module
