module GenericRefFEs

using Gridap
using Gridap.Helpers

using LinearAlgebra

export GenericRefFE
export RTRefFE
export NedelecRefFE
export GenericDOFBasis

import Gridap: shfbasis
import Gridap: polytope
import Gridap: nfacedofs
import Gridap: dofbasis
import Gridap: evaluate!

import Base:length

# """
# It represents the DOF basis for a generic reference FE, in which we provide the vector of nodes
# in which a function must be evaluated in order to be able to compute a set of
# moments. The array of moments (functionals) take the vector of nodal values of
# the function and returns the evaluation of moments. Thus, we can define the
# set of nodes as $\{ xn_i \}$, and the set of moments as $m_i = \sum m_ij f(xn_j)$.
# We note that $m_ij$ is a field value, i.e., if $f$ is a scalar (resp., vector,
# or tensor) field then $m_ij$ is a scalar, (resp., vector, or tensor) value.
#
# We note that the nodes and moments arrays are stored n-face wise. Nodes and
# moments in reference FEs are owned by a particular n-face. By imposing continuity
# of these moments (degrees of freedom) across cells in the physical space we can
# prove a particular type of continuity (e.g., normal, tangent, or full continuity
# for div, curl, or grad-conforming FE spaces). Since the nodal values required
# for every n-face are different, to segregate these computations n-face-wise does
# not introduce any overhead. In fact, it is more effective for the storage of
# $m_ij$, since one does not have to store the zeros related to nodes of other
# n-faces.
# """
struct GenericDOFBasis{D,T,S} <: DOFBasis{D,T}
  nodes::Vector{Array{Point{D,S}}}
  moments::Vector{Array{T}}
  _cache_field#::Vector{T}
  _cache_basis#::Matrix{T}
end

nodes(b::GenericDOFBasis) = b.nodes

moments(b::GenericDOFBasis) = b.moments

length(b::GenericDOFBasis{D,T} where {D,T}) = sum([size(i,1) for i in b.moments])

# """
# It represents the reference FE related to the `GenericDOFBasis` (see
# documentation for more details). The rest of ingredients are standard. It
# involves a polytope, a basis generated as a prebasis and change-of-basis to
# get the canonical basis related to the DOF basis, and an array of nfacedofs
# that provide the ids of dofs for a given n-face.
# """
# @santiagobadia : Do we really need nfacedofs?
# @santiagobadia : We could probably put polytope reference in dofbasis since
# the dofbasis is n-face wise
struct GenericRefFE{D,T} <: RefFE{D,T}
  polytope::Polytope{D}
  dof_basis::GenericDOFBasis{D,T}
  shfbasis::Basis{D,T}
  nfacedofs::Vector{Vector{Int}}
end

dofbasis(this::GenericRefFE{D,T} where {D,T})::DOFBasis = this.dof_basis

polytope(this::GenericRefFE{D,T} where {D,T})::Polytope = this.polytope

shfbasis(this::GenericRefFE{D,T} where {D,T})::Basis = this.shfbasis

nfacedofs(this::GenericRefFE{D,T} where {D,T})::Vector{Vector{Int}} = this.nfacedofs

function _GenericRefFE(p::Polytope{D}, dof_basis::GenericDOFBasis{D,T,S}, shfbasis::Basis{D,T}, nfacedofs) where {D,T,S}
  GenericRefFE{D,T}(p,dof_basis,shfbasis,nfacedofs)
end

function GenericDOFBasis(nodes::Vector{Array{Point{D,S}}}, moments::Vector{Array{T}}) where {D,T,S}
  ndofs = sum([size(i,1) for i in moments])
  nnodes = [length(i) for i in nodes]
  cache_field = [ zeros(T,i) for i in nnodes]
  cache_basis = [ zeros(T,ndofs,i) for i in nnodes]

  GenericDOFBasis{D,T,S}(
  nodes,
  moments,
  cache_field,
  cache_basis)
end

function evaluate!(
  b::GenericDOFBasis{D,T},f::Field{D,T},dofs::AbstractVector{E}) where {D,T,E}
  for (n,v) in zip(b.nodes,b._cache_field)
    evaluate!(f,n,v)
  end
  k = 0
  for (m,v) in zip(b.moments,b._cache_field)
    l = size(m,1)
    if length(m) > 0
      dofs[k+1:k+l] = m*v
      k += l
    end
  end
  return dofs
end

function evaluate!(
  b::GenericDOFBasis{D,T},f::Basis{D,T},dofs::AbstractMatrix{E}) where {D,T,E}
  for (n,v) in zip(b.nodes,b._cache_basis)
    evaluate!(f,n,v)
  end
  k = 0
  for (m,v) in zip(b.moments,b._cache_basis)
    l = size(m,1)
    if length(m) > 0
      dofs[:,k+1:k+l] = v*m'
      k += l
    end
  end
  dofs
end

function _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = ConstantCellValue(fips,nc)
  c_wips = ConstantCellValue(wips,nc)
  pquad = evaluate(fgeomap,c_fips)
  c_fips, pquad, c_wips
end

# Ref FE to faces geomaps
function _ref_face_to_faces_geomap(p,fp)
  cfvs = nfaces_vertices(p,dim(fp))
  nc = length(cfvs)
  freffe = LagrangianRefFE(Float64,fp,1)
  fshfs = shfbasis(freffe)
  cfshfs = ConstantCellValue(fshfs, nc)
  fgeomap = lincomb(cfshfs,cfvs)
end

function _broadcast(::Type{T},n,b) where T
  c = Array{T}(undef,size(b))
  for (ii, i) in enumerate(b)
    c[ii] = i*n
  end
  return c
end

function _broadcast_cross(::Type{T},n,b) where T
  c = Array{T}(undef,size(b))
  for (ii, i) in enumerate(b)
    c[ii] = T(cross(i.array,n.array))# cross product
  end
  return c
end

function _broadcast_extend(::Type{T},Tm,b) where T
  c = Array{T}(undef,size(b))
  for (ii,i) in enumerate(b)
    c[ii] = T(Tm*[i...])
  end
  return c
end

function _nfacedofs_basis(p,moments)
  ndofs = [size(moments[i],1) for i in 1:length(moments)]
  _nfacedofs = Vector{Vector{Int}}(undef,num_nfaces(p))
  _nfacedofs[1] = Int[1:ndofs[1]...]
  k = ndofs[1]
  for i in 2:length(ndofs)
    _nfacedofs[i] = [k+1:k+ndofs[i]...]
    k += ndofs[i]
  end
  return _nfacedofs
end

function _nfaces_array_dim!(p,dim,array,nf_vals)
  nfs = nfaces_dim(p,dim)
  array[nfs] = nf_vals
end

function _initialize_arrays(prebasis,p)

  # Field, point, and entry types
  ft = eltype(prebasis)
  et = eltype(ft)
  pt = Point{dim(p),Float64}

  # Arrays of moments (as its evaluation for all prebasis shape functions)
  # and evaluation points per n-face
  nf_moments = Vector{Array{ft}}(undef,num_nfaces(p))
  nf_nodes = Vector{Array{pt}}(undef,num_nfaces(p))
  nshfs = length(prebasis)
  pb_moments = zeros(et,nshfs,0)

  # Initialize to zero arrays
  zero_moments = zeros(ft,0,0)
  zero_nodes = zeros(pt,0)
  for dim in 0:dim(p)
    for inf in nfaces_dim(p,dim)
      nf_moments[inf] = zero_moments
      nf_nodes[inf] = zero_nodes
    end
  end
  return nf_nodes, nf_moments, pb_moments

end

function _GenericRefFE(p::Polytope,prebasis::Basis,
                      nf_nodes,nf_moments,pb_moments)

  # Change of basis matrix, inv([DF,DC])
  cob = inv(hcat(pb_moments))
  basis = change_basis(prebasis,cob)

  nfacedofs = _nfacedofs_basis(p,nf_moments)

  # Build DOFBasis and RefFE with all this
  dof_basis = GenericDOFBasis(nf_nodes, nf_moments)

  divreffe = _GenericRefFE(p,dof_basis,basis,nfacedofs)

end

function _insert_nface_values!(nf_nodes, nf_moments, pb_moments,
                               prebasis, fcips, fmoments, p, dim)

  # Evaluate basis in faces points, i.e., S(Fi)_{ab} = ϕ^a(xgp_Fi^b)
  pbasis_fcips = [evaluate(prebasis,ps) for ps in fcips]

  # Face moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
  fms_preb = [bps*ms' for (bps,ms) in zip(pbasis_fcips,fmoments)]

  _nfaces_array_dim!(p,dim,nf_moments,fmoments)
  _nfaces_array_dim!(p,dim,nf_nodes,fcips)
  pb_moments = hcat(pb_moments,fms_preb...)

  return nf_nodes, nf_moments, pb_moments

end

include("DivRefFEs.jl")

include("CurlRefFEs.jl")

end # module
