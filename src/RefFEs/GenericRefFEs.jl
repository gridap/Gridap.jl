module GenericRefFEs

using Gridap
using Gridap.Helpers

using LinearAlgebra

export GenericRefFE
export RTRefFE
export NedelecRefFE
export GenericDOFBasis

export field_type

import Gridap: shfbasis
import Gridap: polytope
import Gridap: nfacedofs
import Gridap: dofbasis
import Gridap: evaluate!

import Base:length

struct GenericDOFBasis{D,T,S} <: DOFBasis{D,T}
  nodes::Vector{Array{Point{D,S}}}
  moments::Vector{Array{T}}
  _cache_field#::Vector{T}
  _cache_basis#::Matrix{T}
end

length(b::GenericDOFBasis{D,T} where {D,T}) = sum([size(i,1) for i in b.moments])

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

function field_type(this::GenericRefFE{D,T}) where {D,T}
  return T
end

function RTRefFE(p:: Polytope, order::Int)

  if !(all(extrusion(p).array .== HEX_AXIS))
    @notimplemented
  end

  # Prebasis
  et = Float64
  prebasis = CurlGradMonomialBasis(VectorValue{dim(p),et},order)


  # Field, point, and entry types
  ft = VectorValue{dim(p),Float64}
  pt = Point{dim(p),Float64}

  # Arrays of moments (as its evaluation for all prebasis shape functions)
  # and evaluation points per n-face
  nface_moments = Vector{Array{ft}}(undef,num_nfaces(p))
  nface_evaluation_points = Vector{Array{pt}}(undef,num_nfaces(p))
  nshfs = length(prebasis)
  preb_eval = zeros(et,nshfs,0)

  if (order == 1)
    dims = [collect(0:dim(p)-2)...,dim(p)]
  else
    dims = [collect(0:dim(p)-2)...]
  end

  _null_nface_dim!(p,dims,et,nface_moments,nface_evaluation_points)

  # Face moments

  # Reference facet
  fp = ref_nface_polytope(p,dim(p)-1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial faces
  degree = order*2
  fquad = Quadrature(fp,degree)
  fips = coordinates(fquad)
  wips = weights(fquad)
  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fshfs = Gridap.RefFEs._monomial_basis(fp,Float64,order-1)
  fmoments = _RT_face_moments(p, fshfs, c_fips, fcips, fwips)

  # Evaluate basis in faces points, i.e., S(Fi)_{ab} = ϕ^a(xgp_Fi^b)
  pbasis_fcips = [evaluate(prebasis,ps) for ps in fcips]

  # Face moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
  fms_preb = [bps*ms' for (bps,ms) in zip(pbasis_fcips,fmoments)]

  _nfaces_array_dim!(p,dim(p)-1,nface_moments,fmoments)
  _nfaces_array_dim!(p,dim(p)-1,nface_evaluation_points,fcips)
  preb_eval = hcat(preb_eval,fms_preb...)

  # Cell moments

  if (order > 1)

    # Compute integration points at interior
    degree = 2*order
    iquad = Quadrature(p,degree)
    ccips = coordinates(iquad)
    cwips = weights(iquad)

    # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
    cbasis = GradMonomialBasis(VectorValue{dim(p),Float64},order-1)
    cmoments = _RT_cell_moments(p, cbasis, ccips, cwips )

    # Evaluate basis in cell points, i.e., S(C)_{ab} = ϕ^a(xgp_C^b)
    pbasis_ccips = evaluate(prebasis,ccips)

    # Cell moments evaluated for basis, i.e., DC = S(C)*M(C)^T
    cms_preb = pbasis_ccips*cmoments'

    _nfaces_array_dim!(p,dim(p),nface_moments,[cmoments])
    _nfaces_array_dim!(p,dim(p),nface_evaluation_points,[ccips])
    preb_eval = hcat(preb_eval,cms_preb)

  end

  # Change of basis matrix, inv([DF,DC])
  cob = inv(hcat(preb_eval))
  basis = change_basis(prebasis,cob)

  nfacedofs = _nfacedofs_basis(p,nface_moments)

  # Build DOFBasis and RefFE with all this
  dof_basis = GenericDOFBasis(nface_evaluation_points, nface_moments)

  divreffe = _GenericRefFE(p,dof_basis,basis,nfacedofs)
end

function NedelecRefFE(p:: Polytope, order::Int)

  if !(all(extrusion(p).array .== HEX_AXIS))
    @notimplemented
  end

  # Prebasis
  et = Float64
  prebasis = GradMonomialBasis(VectorValue{dim(p),et},order)

  # Field, point, and entry types
  ft = VectorValue{dim(p),Float64}
  pt = Point{dim(p),Float64}

  # Arrays of moments (as its evaluation for all prebasis shape functions)
  # and evaluation points per n-face
  nface_moments = Vector{Array{ft}}(undef,num_nfaces(p))
  nface_evaluation_points = Vector{Array{pt}}(undef,num_nfaces(p))
  nshfs = length(prebasis)
  preb_eval = zeros(et,nshfs,0)

  if (order == 1)
    dims = [0,collect(2:dim(p))...]
  else
    dims = [0]
  end

  println(dims)
  _null_nface_dim!(p,dims,et,nface_moments,nface_evaluation_points)

  # Edge moments (dim = 1)

  # Reference facet
  dim1 = 1
  fp = ref_nface_polytope(p,dim1) # dim 1

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial edges
  degree = order*2
  fquad = Quadrature(fp,degree)
  fips = coordinates(fquad)
  wips = weights(fquad)
  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Edge moments, i.e., M(Ei)_{ab} = q_RE^a(xgp_REi^b) w_Fi^b t_Ei ⋅ ()
  fshfs = Gridap.RefFEs._monomial_basis(fp,Float64,order-1)
  fmoments = _Nedelec_edge_moments(p, fshfs, c_fips, fcips, fwips)

  # Evaluate basis in edges points, i.e., S(Ei)_{ab} = ϕ^a(xgp_Ei^b)
  pbasis_fcips = [evaluate(prebasis,ps) for ps in fcips]

  # Edge moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
  fms_preb = [bps*ms' for (bps,ms) in zip(pbasis_fcips,fmoments)]

  _nfaces_array_dim!(p,1,nface_moments,fmoments)
  _nfaces_array_dim!(p,1,nface_evaluation_points,fcips)
  preb_eval = hcat(preb_eval,fms_preb...)

  # Nedelec Face moments (dim = 1)

  if ( dim(p) == 3 && order > 1)
    # Reference facet
    dimf = dim(p)-1
    fp = ref_nface_polytope(p,dimf) # dim 1

    # geomap from ref face to polytope faces
    fgeomap = _ref_face_to_faces_geomap(p,fp)

    # Compute integration points at all polynomial edges
    degree = order*2
    fquad = Quadrature(fp,degree)
    fips = coordinates(fquad)
    wips = weights(fquad)
    c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

    # Face moments, i.e., M(Fi)_{ab} = w_Fi^b q_RF^a(xgp_RFi^b) (n_Fi × ())
    fshfs = GradMonomialBasis(VectorValue{dim(fp),et},order-1)
    fmoments = _Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)

    # Evaluate basis in edges points, i.e., S(Ei)_{ab} = ϕ^a(xgp_Ei^b)
    pbasis_fcips = [evaluate(prebasis,ps) for ps in fcips]

    # Edge moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
    fms_preb = [bps*ms' for (bps,ms) in zip(pbasis_fcips,fmoments)]

    _nfaces_array_dim!(p,dim(p)-1,nface_moments,fmoments)
    _nfaces_array_dim!(p,dim(p)-1,nface_evaluation_points,fcips)
    preb_eval = hcat(preb_eval,fms_preb...)
  end
  # Cell moments

  if (order > 1)

    # Compute integration points at interior
    degree = 2*order
    iquad = Quadrature(p,degree)
    ccips = coordinates(iquad)
    cwips = weights(iquad)

    # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
    cbasis = CurlGradMonomialBasis(VectorValue{dim(p),Float64},order-1)
    cmoments = _Nedelec_cell_moments(p, cbasis, ccips, cwips )

    # Evaluate basis in cell points, i.e., S(C)_{ab} = ϕ^a(xgp_C^b)
    pbasis_ccips = evaluate(prebasis,ccips)

    # Cell moments evaluated for basis, i.e., DC = S(C)*M(C)^T
    cms_preb = pbasis_ccips*cmoments'

    _nfaces_array_dim!(p,dim(p),nface_moments,[cmoments])
    _nfaces_array_dim!(p,dim(p),nface_evaluation_points,[ccips])
    preb_eval = hcat(preb_eval,cms_preb)

  end

  # Change of basis matrix, inv([DF,DC])
  cob = inv(hcat(preb_eval))
  basis = change_basis(prebasis,cob)

  nfacedofs = _nfacedofs_basis(p,nface_moments)

  # Build DOFBasis and RefFE with all this
  dof_basis = GenericDOFBasis(nface_evaluation_points, nface_moments)

  divreffe = _GenericRefFE(p,dof_basis,basis,nfacedofs)
end

function _GenericRefFE(p::Polytope{D}, dof_basis::GenericDOFBasis{D,T,S}, shfbasis::Basis{D,T}, nfacedofs) where {D,T,S}
  GenericRefFE{D,T}(p,dof_basis,shfbasis,nfacedofs)
end

# function GenericDOFBasis(nodes::Vector{Array{Point{D,T}}}, moments::Vector{Array{Point{D,T}}}) where {D,T}
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

function _RT_cell_moments(p, cbasis, ccips, cwips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips'.*ishfs_iips
end

const _Nedelec_cell_moments = _RT_cell_moments

# Ref facet FE functions evaluated at the facet integration points (in ref facet)
function _RT_face_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  cfshfs = ConstantCellValue(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  fns, os = face_normals(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  cvals = [ _broadcast(typeof(n),n*o,b) for (n,o,b) in zip(fns,os,cvals)]
  return cvals
end

# Ref facet FE functions evaluated at the facet integration points (in ref facet)
function _Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)
  nc = length(c_fips)
  fvs = nfaces_vertices(p,dim(p)-1)
  fts = [hcat([vs[2]-vs[1]...],[vs[3]-vs[1]...]) for vs in fvs]
  cfshfs = ConstantCellValue(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  fns, os = face_normals(p)
  # @santiagobadia : Temporary hack for making it work for structured hex meshes
  ft = eltype(fns)
  cvals = [ Gridap.RefFEs.GenericRefFEs._broadcast_extend(ft,Tm,b) for (Tm,b) in zip(fts,cvals)]
  cvals = [ Gridap.RefFEs.GenericRefFEs._broadcast_cross(ft,n,b) for (n,b) in zip(fns,cvals)]
  return cvals
end

function _Nedelec_edge_moments(p, fshfs, c_fips, fcips, fwips)
  ts = edge_tangents(p)
  nc = length(c_fips)
  cfshfs = ConstantCellValue(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  # @santiagobadia : Only working for oriented meshes now
  cvals = [ _broadcast(typeof(t),t,b) for (t,b) in zip(ts,cvals)]
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
function _ref_face_to_faces_geomap(p,fp)
  cfvs = nfaces_vertices(p,dim(fp))
  nc = length(cfvs)
  freffe = LagrangianRefFE(Float64,fp,1)
  fshfs = shfbasis(freffe)
  cfshfs = ConstantCellValue(fshfs, nc)
  fgeomap = lincomb(cfshfs,cfvs)
end

# Ref FE to faces geomaps
function _geomap_jacobian(p,fp)
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

function _nfaces_array_dim!(p,dim,array,nf_vals)
  nfs = nfaces_dim(p,dim)
  array[nfs] = nf_vals
end

function _null_nface_dim!(p,dims,et,nface_moments,nface_evaluation_points)

  ft = VectorValue{dim(p),Float64}
  pt = Point{dim(p),Float64}
  zero_moments = zeros(ft,0,0)
  zero_ips = zeros(pt,0)
  for dim in dims
    for inf in nfaces_dim(p,dim)
      nface_moments[inf] = zero_moments
      nface_evaluation_points[inf] = zero_ips
    end
  end

end

end # module
