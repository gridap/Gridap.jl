
##
using Gridap

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

function _nfaces_array(fmoments,cmoments)
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
  # Facet ref polytope
################################################################################
p = Polytope(1,1)
order = 2

fp = _ref_face_polytope(p)

c_fvs = _face_vertices(p)

fgeomap = _ref_face_to_faces_geomap(p,fp,c_fvs)

# Compute integration points at all polynomial faces
c_fips, fcips, fwips = _faces_quadrature_points_weights(p, fp, fgeomap, order)

# Compute integration points at interior
ccips, cwips = _cell_quadrature_points_weights(p, order)

# Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
fmoments = _face_moments(p, fp, order, c_fips, fwips, c_fvs )

# Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
if(order>1) cmoments = _cell_moments(p, order, ccips, cwips ) end

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
if (order > 1) cms_preb = pbasis_ccips*cmoments' end

# Change of basis matrix, inv([DF,DC])
order > 1 ? cob = inv(hcat(fms_preb...,cms_preb)) : cob = inv(hcat(fms_preb...))

# Store S(NF) and M(NF) for all NFs of polytope
SNF = _nfaces_array(fmoments,cmoments)

# I want to store [Xgp_NF] for all NFS of polytope "Nodes per n-face"
IPNF = _nfaces_array(fcips,ccips)
isa(SNF,Vector{Array{Point{dim(p),Float64}}})
typeof(SNF)

# Compute nfacedofs (prev version)
# Build DOFBasis and RefFE with all this

##

struct DivConformingRefFE{D,T} <: RefFE{D,T}
  polytope::Polytope{D}
  dof_basis::DivComformingDOFBasis{D,T}
  # dof basis TO BE DONE, everything is already in the implementation...
  # Instead of evaluating in nodes, we should evaluate somewhere else,
  # i.e., the facet/interior integration points
  shfbasis::Basis{D,T}
  nfacedofs::Vector{Vector{Int}}
end

struct DivConformingDOFBasis{D,T} <: DOFBasis{D,T}
  nodes::Vector{Array{Point{dim(p),T}}})
  moments::Vector{Array{Point{dim(p),T}}})
  # dof_to_node::Vector{Int} # Non-sense anymore
  # dof_to_comp::Vector{Int} # Non-sense anymore
  # node_and_comp_to_dof::Matrix{Int} # Non-sense anymore
  _cache_field::Vector{T}
  _cache_basis::Matrix{T}
  # Missing part weights
end
##
