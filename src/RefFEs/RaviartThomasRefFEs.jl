module RaviartThomasRefFEs

using Gridap
using Gridap.Helpers

export RaviartThomasRefFE

import Gridap: shfbasis
import Gridap: polytope
import Gridap: nfacedofs
import Gridap: dofbasis

struct RaviartThomasDOFBasis{D,T,S} <: DOFBasis{D,T}
  nodes::Vector{Array{Point{D,S}}}
  moments::Vector{Array{T}}
  _cache_field::Vector{S}
  _cache_basis::Matrix{S}
end

struct RaviartThomasRefFE{D,T} <: RefFE{D,T}
  polytope::Polytope{D}
  dof_basis::RaviartThomasDOFBasis{D,T}
  shfbasis::Basis{D,T}
  nfacedofs::Vector{Vector{Int}}
end

dofbasis(this::RaviartThomasRefFE{D,T} where {D,T})::DOFBasis = this.dof_basis

polytope(this::RaviartThomasRefFE{D,T} where {D,T})::Polytope = this.polytope

shfbasis(this::RaviartThomasRefFE{D,T} where {D,T})::Basis = this.shfbasis

nfacedofs(this::RaviartThomasRefFE{D,T} where {D,T})::Vector{Vector{Int}} = this.nfacedofs

function RaviartThomasRefFE(p:: Polytope, order::Int)

  if !(all(extrusion(p).array .== HEX_AXIS))
    @notimplemented
  end

  # Reference facet
  fp = _ref_face_polytope(p)

  c_fvs = _face_vertices(p)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp,c_fvs)

  # Compute integration points at all polynomial faces
  c_fips, fcips, fwips = _faces_quadrature_points_weights(p, fp, fgeomap, order)

  # Compute integration points at interior
  ccips, cwips = _cell_quadrature_points_weights(p, order)

  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fmoments = _face_moments(p, fp, order, c_fips, fwips, c_fvs )

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  if (order > 1)
    cmoments = _cell_moments(p, order, ccips, cwips )
  else
    cmoments = eltype(eltype(fmoments))[]
  end

  # Prebasis
  prebasis = CurlGradMonomialBasis(VectorValue{dim(p),Float64},order)

  # Evaluate basis in faces points, i.e., S(Fi)_{ab} = ϕ^a(xgp_Fi^b)
  nc = num_nfaces(p,dim(p)-1)
  c_prebasis = ConstantCellValue(prebasis, nc)
  pbasis_fcips = evaluate(c_prebasis,fcips)

  # Evaluate basis in cell points, i.e., S(C)_{ab} = ϕ^a(xgp_C^b)
  if(order>1) pbasis_ccips = evaluate(prebasis,ccips) end

  # Face moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
  fms_preb = [pbasis_fcips[i]*fmoments[i]' for i in 1:nc]

  # Cell moments evaluated for basis, i.e., DC = S(C)*M(C)^T
  if (order > 1)
    cms_preb = pbasis_ccips*cmoments'
  else
    cms_preb = eltype(eltype(fms_preb))[]
  end

  # Change of basis matrix, inv([DF,DC])
  order > 1 ? cob = inv(hcat(fms_preb...,cms_preb)) : cob = inv(hcat(fms_preb...))
  basis = change_basis(prebasis,cob)

  # Store S(NF) and M(NF) for all NFs of polytope
  moments = _nfaces_array(p,fmoments,cmoments)

  # I want to store [Xgp_NF] for all NFS of polytope "Nodes per n-face"
  int_points = _nfaces_array(p,fcips,ccips)

  nfacedofs = _nfacedofs_basis(p,moments)

  # Build DOFBasis and RefFE with all this
  dof_basis = RaviartThomasDOFBasis(int_points, moments)

  divreffe = _RaviartThomasRefFE(p,dof_basis,basis,nfacedofs)
end

function _RaviartThomasRefFE(p::Polytope{D}, dof_basis::RaviartThomasDOFBasis{D,T,S}, shfbasis::Basis{D,T}, nfacedofs) where {D,T,S}
  RaviartThomasRefFE{D,T}(p,dof_basis,shfbasis,nfacedofs)
end

# function RaviartThomasDOFBasis(nodes::Vector{Array{Point{D,T}}}, moments::Vector{Array{Point{D,T}}}) where {D,T}
function RaviartThomasDOFBasis(nodes::Vector{Array{Point{D,S}}}, moments::Vector{Array{T}}) where {D,T,S}
  ndofs = sum(size(moments[i],1) for i in 1:length(moments))
  cache_field = zeros(S,ndofs)
  cache_basis = zeros(S,ndofs,ndofs) # TO CHECK !!!

  RaviartThomasDOFBasis{D,T,S}(
  nodes,
  moments,
  cache_field,
  cache_basis)
end


function _cell_moments(p, order, ccips, cwips)
  cbasis = GradMonomialBasis(VectorValue{dim(p),Float64},order-1)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(cbasis,ccips)
  return cwips'.*ishfs_iips
end

# Ref facet FE functions evaluated at the facet integration points (in ref facet)
function _face_moments(p, fp, order, c_fips, fwips, c_fvs)
  nc = num_nfaces(p,dim(p)-1)
  fshfs = Gridap.RefFEs._monomial_basis(fp,Float64,order-1)
  cfshfs = ConstantCellValue(fshfs, nc)
  cvals = evaluate(cfshfs,c_fips)
  cvals = [fwips[i]'.*cvals[i] for i in 1:nc]
  fns, o = facet_normals(p)
  cvals = [ _mybroadcast(typeof(n),n,b) for (n,b) in zip(fns,cvals)]
  return cvals
end

function _cell_quadrature_points_weights(p, order)
  iquad = Quadrature(p,order*2)
  iips = coordinates(iquad)
  iwei = weights(iquad)
  iips, iwei
end

function _faces_quadrature_points_weights(p, fp, fgeomap, order)
  nc = num_nfaces(p,dim(p)-1)
  fquad = Quadrature(fp,order*2)
  fips = coordinates(fquad)
  c_fips = ConstantCellValue(fips,nc)
  cquad = ConstantCellQuadrature(fquad,nc)
  xquad = coordinates(cquad)
  wquad = weights(cquad)
  pquad = evaluate(fgeomap,xquad)
  c_fips, pquad, wquad
end

# Ref FE to faces geomaps
function _ref_face_to_faces_geomap(p,fp,cfvs)
  nc = num_nfaces(p,dim(p)-1)
  freffe = LagrangianRefFE(Float64,fp,1)
  fshfs = shfbasis(freffe)
  cfshfs = ConstantCellValue(fshfs, nc)
  fgeomap = lincomb(cfshfs,cfvs)
end

function _face_vertices(p)
  nc = num_nfaces(p,dim(p)-1)
  verts = vertices_coordinates(p)
  faces_vs = nface_connections(p,dim(p)-1,0)
  fvs = Gridap.CellValuesGallery.CellValueFromArray(faces_vs)
  vs = Gridap.CellValuesGallery.CellValueFromArray(verts)
  cfvs = Gridap.CellValuesGallery.CellVectorFromLocalToGlobal(fvs,vs)
end

function _mybroadcast(::Type{T},n,b) where T
  c = Array{T}(undef,size(b))
  for (ii, i) in enumerate(b)
    c[ii] = i*n
  end
  return c
end

# Check all facets are of the same typefps = nface_ref_polytopes(p)[facets]
function _ref_face_polytope(p)
  faces = nfaces_dim(p,dim(p)-1)
  fps = nface_ref_polytopes(p)[faces]
  @assert(all(extrusion.(fps) .== extrusion(fps[1])), "All facet must be of the same type")
  return fps[1]
end

function _nfaces_array(p,fmoments,cmoments)
  T = eltype(eltype(fmoments))
  SNF = Vector{Array{T}}(undef,num_nfaces(p))
  zeromat = eltype(eltype(fmoments))[]
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

end # module
