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

function NormalVector(
  phi::CellField,
  desc::BoundaryDescriptor,
  polytopes::CellValue{<:Polytope})
  @notimplemented
end

function NormalVector(
  phi::CellField{D},
  desc::BoundaryDescriptor,
  polytopes::ConstantCellValue{<:Polytope{D}}) where D

  jac = gradient(phi)
  jac_surf = restrict(jac,desc)

  p = polytopes.value
  n, _ = facet_normals(p)

  _n = [ [ni,] for ni in n  ]

  facet_to_lfacet = desc.facet_to_lfacet

  nvec_ref = IndexCompressedCellValue(_n,facet_to_lfacet)

  J = typeof(jac_surf)
  R = typeof(nvec_ref)

  Z = D -1
  NormalVector{Z,D,J,R}(jac_surf,nvec_ref)

end

length(nv::NormalVector) = length(jac_surf)

function evaluate(nv::NormalVector{Z}, q::CellPoints{Z}) where Z
  jac_surf_q = evaluate(nv.jac_surf,q)
  apply(_map_normal,jac_surf_q,nv.nvec_ref,broadcast=true)
end

function _map_normal(J,n)
  v = inv(J')*n
  m = sqrt(inner(v,v))
  if m < eps()
    return zero(n)
  else
    return v/m
  end
end

end # module
