module NewDivRefFEs

using Gridap
using Gridap.Helpers

export NewDivRefFE
export NewDivDOFBasis

export field_type

import Gridap: shfbasis
import Gridap: polytope
import Gridap: nfacedofs
import Gridap: dofbasis
import Gridap: evaluate!

import Base:length

struct NewDivDOFBasis{D,T,S} <: DOFBasis{D,T}
  nodes::Vector{Array{Point{D,S}}}
  moments::Vector{Array{T}}
  _cache_field#::Vector{T}
  _cache_basis#::Matrix{T}
end

length(b::NewDivDOFBasis{D,T} where {D,T}) = sum([size(i,1) for i in b.moments])

struct NewDivRefFE{D,T} <: RefFE{D,T}
  polytope::Polytope{D}
  dof_basis::NewDivDOFBasis{D,T}
  shfbasis::Basis{D,T}
  nfacedofs::Vector{Vector{Int}}
end

dofbasis(this::NewDivRefFE{D,T} where {D,T})::DOFBasis = this.dof_basis

polytope(this::NewDivRefFE{D,T} where {D,T})::Polytope = this.polytope

shfbasis(this::NewDivRefFE{D,T} where {D,T})::Basis = this.shfbasis

nfacedofs(this::NewDivRefFE{D,T} where {D,T})::Vector{Vector{Int}} = this.nfacedofs

function field_type(this::NewDivRefFE{D,T}) where {D,T}
  return T
end

function NewDivRefFE(p:: Polytope, order::Int)

  if !(all(extrusion(p).array .== HEX_AXIS))
    @notimplemented
  end

  # Prebasis
  prebasis = CurlGradMonomialBasis(VectorValue{dim(p),Float64},order)
  ft = VectorValue{dim(p),Float64}
  pt = Point{dim(p),Float64}
  et = eltype(ft)

  nshfs = length(prebasis)

  # Face moments

  # Reference facet
  fp = ref_nface_polytope(p,dim(p)-1)
  c_fvs = nfaces_vertices(p,dim(p)-1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp,c_fvs)

  # Compute integration points at all polynomial faces
  degree = order*2
  fquad = Quadrature(fp,degree)
  fips = coordinates(fquad)
  wips = weights(fquad)
  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fshfs = Gridap.RefFEs._monomial_basis(fp,Float64,order-1)
  fmoments = _nface_moments(p, fshfs, c_fips, fcips, fwips)

  # Evaluate basis in faces points, i.e., S(Fi)_{ab} = ϕ^a(xgp_Fi^b)
  nc = length(fcips)
  c_prebasis = ConstantCellValue(prebasis, nc)
  pbasis_fcips = evaluate(c_prebasis,fcips)

  # Face moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
  fms_preb = [pbasis_fcips[i]*fmoments[i]' for i in 1:nc]

  # Cell moments

  if (order > 1)

    # Compute integration points at interior
    degree = 2*order
    iquad = Quadrature(p,degree)
    ccips = coordinates(iquad)
    cwips = weights(iquad)

    # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
    cbasis = GradMonomialBasis(VectorValue{dim(p),Float64},order-1)
    cmoments = _cell_moments(p, cbasis, ccips, cwips )

    # Evaluate basis in cell points, i.e., S(C)_{ab} = ϕ^a(xgp_C^b)
    pbasis_ccips = evaluate(prebasis,ccips)

    # Cell moments evaluated for basis, i.e., DC = S(C)*M(C)^T
    cms_preb = pbasis_ccips*cmoments'

  else

    cmoments = zeros(ft,0,0)

    cms_preb = zeros(et,nshfs,0)

    ccips = zeros(pt,0)

  end

  # Change of basis matrix, inv([DF,DC])
  # order > 1 ? cob = inv(hcat(fms_preb...,cms_preb)) : cob = inv(hcat(fms_preb...))
  cob = inv(hcat(fms_preb...,cms_preb))
  basis = change_basis(prebasis,cob)

  # Store S(NF) and M(NF) for all NFs of polytope
  moments = _nfaces_array(p,fmoments,cmoments,ft)

  # I want to store [Xgp_NF] for all NFS of polytope "Nodes per n-face"
  int_points = _nfaces_array(p,fcips,ccips,pt)

  nfacedofs = _nfacedofs_basis(p,moments)

  # Build DOFBasis and RefFE with all this
  dof_basis = NewDivDOFBasis(int_points, moments)

  divreffe = _NewDivRefFE(p,dof_basis,basis,nfacedofs)
end

function _NewDivRefFE(p::Polytope{D}, dof_basis::NewDivDOFBasis{D,T,S}, shfbasis::Basis{D,T}, nfacedofs) where {D,T,S}
  NewDivRefFE{D,T}(p,dof_basis,shfbasis,nfacedofs)
end

# function NewDivDOFBasis(nodes::Vector{Array{Point{D,T}}}, moments::Vector{Array{Point{D,T}}}) where {D,T}
function NewDivDOFBasis(nodes::Vector{Array{Point{D,S}}}, moments::Vector{Array{T}}) where {D,T,S}
  ndofs = sum([size(i,1) for i in moments])
  nnodes = [length(i) for i in nodes]
  cache_field = [ zeros(T,i) for i in nnodes]
  cache_basis = [ zeros(T,ndofs,i) for i in nnodes]

  NewDivDOFBasis{D,T,S}(
  nodes,
  moments,
  cache_field,
  cache_basis)
end

function evaluate!(
  b::NewDivDOFBasis{D,T},f::Field{D,T},dofs::AbstractVector{E}) where {D,T,E}
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
  b::NewDivDOFBasis{D,T},f::Basis{D,T},dofs::AbstractMatrix{E}) where {D,T,E}
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

function _cell_moments(p, cbasis, ccips, cwips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips'.*ishfs_iips
end

# Ref facet FE functions evaluated at the facet integration points (in ref facet)
function _nface_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  cfshfs = ConstantCellValue(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  fns, os = facet_normals(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  cvals = [ _broadcast(typeof(n),n*o,b) for (n,o,b) in zip(fns,os,cvals)]
  return cvals
end

function _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)
  nc = length(fgeomap)
  c_fips = ConstantCellValue(fips,nc)
  c_wips = ConstantCellValue(wips,nc)
  pquad = evaluate(fgeomap,c_fips)
  c_fips, pquad, c_wips
end

# Ref FE to faces geomaps
function _ref_face_to_faces_geomap(p,fp,cfvs)
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


function _nfaces_array(p,fmoments,cmoments,T)
  SNF = Vector{Array{T}}(undef,num_nfaces(p))
  zeromat = T[]
  for idim in 0:dim(p)-2
    for inf in nfaces_dim(p,idim)
      SNF[inf] = zeromat
    end
  end
  faces = nfaces_dim(p,dim(p)-1)
  SNF[faces] = fmoments
  SNF[end] = cmoments
  return SNF
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

function _null_nface_dim!(p,dim,array,zeros)
  for inf in nfaces_dim(p,dim)
    array[inf] = zeros
  end
end

function _nfaces_array_dim!(p,dim,array,nf_vals)
  nfs = nfaces_dim(p,dim)
  array[nfs] = nf_vals
end

end # module
