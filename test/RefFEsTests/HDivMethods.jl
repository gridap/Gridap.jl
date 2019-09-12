
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

function _ref_facet_polytope(p)
################################################################################
p = Polytope(1,1)
order = 1

# Check all facets are of the same typefps = nface_ref_polytopes(p)[facets]
fp = _ref_facet_polytope(p)

faces = nfaces_dim(p,dim(p)-1)

fps = nface_ref_polytopes(p)[faces]
@assert(all(extrusion.(fps) .== extrusion(fps[1])), "All facet must be of the same type")
c_fvs = _face_vertices(p)
# Facet ref polytope
fp = fps[1]
fgeomap = _ref_face_to_faces_geomap(p,fp,c_fvs)
# Compute integration points at all polynomial faces
c_fips, fcips, fwips = _faces_quadrature_points_weights(p, fp, fgeomap, order)
# Compute integration points at interior
ccips, cwips = _cell_quadrature_points_weights(p, order)
# Face moments
fmoments = _face_moments(p, fp, order, c_fips, fwips, c_fvs )
# b
nfmoments = _face_moments(p, fp, order, c_fips, fwips, c_fvs )
# Cell moments
if(order>1) cmoments = _cell_moments(p, order, ccips, cwips ) end
# Prebasis
prebasis = CurlGradMonomialBasis(VectorValue{dim(p),Float64},order)
# Evaluate basis in faces points
nc = num_nfaces(p,dim(p)-1)
c_prebasis = ConstantCellValue(prebasis, nc)
pbasis_fcips = evaluate(c_prebasis,fcips)
if(order>1) pbasis_ccips = evaluate(prebasis,ccips) end
# Face moments evaluated for basis
fms_preb = [pbasis_fcips[i]*fmoments[i]' for i in 1:nc]
# Cell moments evaluated for basis
if (order > 1)
  cms_preb = pbasis_ccips*cmoments'
  cob = inv(hcat(fms_preb...,cms_preb))
else
  cob = inv(hcat(fms_preb...))
end
# Evaluate basis in cell points
##
