module PDiscRefFEs

using Gridap

export PDiscRefFE

import Gridap: shfbasis
import Gridap: polytope
import Gridap: nfacedofs
import Gridap: dofbasis

struct PDiscRefFE{D,T} <: RefFE{D,T}
  reffe_tet::LagrangianRefFE{D,T}
  polytope::Polytope{D}
  nfacedofs::Vector{Vector{Int}}
end

function PDiscRefFE(::Type{T},D::Integer,order::Integer) where T

  extrusion = fill(HEX_AXIS,D)
  polytope = Polytope(extrusion...)

  PDiscRefFE(T,polytope,order)

end

function PDiscRefFE(
  ::Type{T}, polytope::Polytope{D}, order::Integer) where {T,D}

  extrusion_tet = fill(TET_AXIS,D)
  polytope_tet = Polytope(extrusion_tet...)

  reffe_tet = LagrangianRefFE(T,polytope_tet,order)

  # We assign all DOFs to the interior of the cell.
  # It is not possible to construct continuous spaces
  # with this RefFE
  nnfaces = length(polytope.nfaces)
  ndofs = length(reffe_tet.shfbasis)
  nfacedofs = [Int[] for i in 1:nnfaces]
  nfacedofs[end] = collect(1:ndofs)

  PDiscRefFE{D,T}(reffe_tet,polytope,nfacedofs)

end

dofbasis(this::PDiscRefFE) = this.reffe_tet.dofbasis

polytope(this::PDiscRefFE) = this.polytope

shfbasis(this::PDiscRefFE) = this.reffe_tet.shfbasis

nfacedofs(this::PDiscRefFE) = this.nfacedofs

end # module
