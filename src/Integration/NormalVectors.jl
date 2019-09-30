module NormalVectors

using Gridap
using Gridap.Helpers
using Gridap.CellFields: NonIterableCellFieldLike

export NormalVector

import Base: length
import Gridap: evaluate

struct NormalVector{Z,D,J,R} <: NonIterableCellFieldLike{Z,VectorValue{D,Float64},1}
  jac_surf::J
  nvec_ref::R
end

function NormalVector(trian::BoundaryTriangulation)
  desc = trian.descriptor
  NormalVector(desc)
end

function NormalVector(trian::SkeletonTriangulation)
  desc = trian.descriptor1
  NormalVector(desc)
end

function NormalVector(desc::BoundaryDescriptor)
  phi = desc.cell_phi
  jac = gradient(phi)
  jac_surf = restrict(jac,desc)
  polytopes = desc.cell_to_polytope
  facet_to_lfacet = desc.facet_to_lfacet
  _normal_vector(jac_surf,polytopes,facet_to_lfacet)
end

length(nv::NormalVector) = length(nv.jac_surf)

function evaluate(nv::NormalVector{Z}, q::CellPoints{Z}) where Z
  jac_surf_q = evaluate(nv.jac_surf,q)
  apply(_map_normal,jac_surf_q,nv.nvec_ref,broadcast=true)
end

function _map_normal(J,n)
  v = inv(J)*n
  m = sqrt(inner(v,v))
  if m < eps()
    return zero(n)
  else
    return v/m
  end
end

function _normal_vector(
  jac_surf,
  polytopes::CellValue{<:Polytope},
  facet_to_lfacet)
  @notimplemented
end

function _normal_vector(
  jac_surf,
  polytopes::ConstantCellValue{<:Polytope{D}},
  facet_to_lfacet) where D

  p = polytopes.value
  n, _ = facet_normals(p)

  _n = [ [ni,] for ni in n  ]

  nvec_ref = IndexCompressedCellValue(_n,facet_to_lfacet)

  J = typeof(jac_surf)
  R = typeof(nvec_ref)

  Z = D -1
  NormalVector{Z,D,J,R}(jac_surf,nvec_ref)

end

#@fverdugo this code works does not work here, but works when placed
# in the tests
import Base: show

function show(io::IO,self::NormalVector{Z,D,J,R}) where {Z,D,J,R}
  print(io,"$(nameof(typeof(self))){$Z,$D} object")
end

function show(io::IO,::MIME"text/plain",self::NormalVector{Z,D,J,R}) where {Z,D,J,R}
  show(io,self)
  print(io,":")
  print(io,"\n physdim: $D")
  print(io,"\n refdim: $Z")
  print(io,"\n ncells: $(length(self))")
end


end # module
