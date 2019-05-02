module Wrappers

using Numa.Polytopes
using Numa.Geometry

import UnstructuredGrids.Core: RefCell

#function RefCell(::Polytope{D}) where D
#
#
#  ndims = D
#  faces::Vector{Vector{Vector{Int}}}
#  facetypes::Vector{Vector{Int}}
#  reffaces::Vector{Vector{RefCell}}
#  coordinates::Array{Float64,2}
#  vtkid::Int
#  vtknodes::Vector{Int}
#
#end

function _faces(polytope)

  nface_to_dim, nface_to_vertices, nface_to_code = (
    _extract_nface_to_info(polytope) )

 dim_to_jface_to_vertices, dim_to_jface_to_code = (
    _sort_nfaces_to_vertices_per_dim(
      nface_to_dim, nface_to_vertices, nface_to_code) )

 (dim_to_jface_to_vertices, dim_to_jface_to_code)
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

end # module Wrappers
