module RaviartThomasRefFEsTests
##

using Gridap
using Gridap.Helpers

# using Gridap.RefFEs.RaviartThomasRefFEs

using LinearAlgebra

using Test

using Gridap.RefFEs.NewDivRefFEs
#####
# Create a general Ref FE constructor, which generalizes the current implementations
# for nodal (Lagrangian) and non-nodal (RaviartThomas) RefFEs.
# * A method that generates the geomap from ref n-face to n-face in polytope
# * A functor to be provided by the user (optionally) for every n-face dimension
#####
##
p = Polytope(1,1)
ref_fe = NewDivRefFE(p,3)
ref_fe2 = RaviartThomasRefFE(p,3)

ref_fe.shfbasis.changeofbasis == ref_fe2.shfbasis.changeofbasis

order = 3

if !(all(extrusion(p).array .== HEX_AXIS))
  @notimplemented
end

# Prebasis
et = Float64
prebasis = CurlGradMonomialBasis(VectorValue{dim(p),et},order)
nshfs = length(prebasis)

# Field, point, and entry types
ft = VectorValue{dim(p),Float64}
pt = Point{dim(p),Float64}

# SNF
nface_moments = Vector{Array{ft}}(undef,num_nfaces(p))
# int_points
nface_evaluation_points = Vector{Array{pt}}(undef,num_nfaces(p))
# moments evaluated for prebasis
preb_eval = zeros(et,nshfs,0)

if (order == 1)
  dims = [collect(0:dim(p)-2)...,dim(p)]
else
  dims = [collect(0:dim(p)-2)...]
end

Gridap.RefFEs.NewDivRefFEs._null_nface_dim!(p,dims,et,nface_moments,nface_evaluation_points)


##
# Faces

# Reference facet
fp = ref_nface_polytope(p,dim(p)-1)
c_fvs = nfaces_vertices(p,dim(p)-1)

# geomap from ref face to polytope faces
fgeomap = Gridap.RefFEs.NewDivRefFEs._ref_face_to_faces_geomap(p,fp,c_fvs)

# Compute integration points at all polynomial faces
degree = order*2
fquad = Quadrature(fp,degree)
fips = coordinates(fquad)
wips = weights(fquad)
c_fips, fcips, fwips = Gridap.RefFEs.NewDivRefFEs._nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

# Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
fshfs = Gridap.RefFEs._monomial_basis(fp,Float64,order-1)
fmoments = Gridap.RefFEs.NewDivRefFEs._nface_moments(p, fshfs, c_fips, fcips, fwips)

# Evaluate basis in faces points, i.e., S(Fi)_{ab} = ϕ^a(xgp_Fi^b)
nc = length(fcips)
c_prebasis = ConstantCellValue(prebasis, nc)
pbasis_fcips = evaluate(c_prebasis,fcips)

# Face moments evaluated for basis, i.e., DF = [S(F1)*M(F1)^T, …, S(Fn)*M(Fn)^T]
fms_preb = [pbasis_fcips[i]*fmoments[i]' for i in 1:nc]

##
Gridap.RefFEs.NewDivRefFEs._nfaces_array_dim!(p,dim(p)-1,nface_moments,fmoments)
Gridap.RefFEs.NewDivRefFEs._nfaces_array_dim!(p,dim(p)-1,nface_evaluation_points,fcips)
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
  cmoments = Gridap.RefFEs.NewDivRefFEs._cell_moments(p, cbasis, ccips, cwips )

  # Evaluate basis in cell points, i.e., S(C)_{ab} = ϕ^a(xgp_C^b)
  pbasis_ccips = evaluate(prebasis,ccips)

  # Cell moments evaluated for basis, i.e., DC = S(C)*M(C)^T
  cms_preb = pbasis_ccips*cmoments'

  Gridap.RefFEs.NewDivRefFEs._nfaces_array_dim!(p,dim(p),nface_moments,[cmoments])
  Gridap.RefFEs.NewDivRefFEs._nfaces_array_dim!(p,dim(p),nface_evaluation_points,[ccips])
  preb_eval = hcat(preb_eval,cms_preb)

else

  Gridap.RefFEs.NewDivRefFEs._null_nface_dim!(p,dim(p),nface_moments,zero_moments)
  Gridap.RefFEs.NewDivRefFEs._null_nface_dim!(p,dim(p),nface_evaluation_points,zero_ips)

end


# Change of basis matrix, inv([DF,DC])
cob = inv(hcat(preb_eval))
basis = change_basis(prebasis,cob)


# Store S(NF) and M(NF) for all NFs of polytope
# moments = Gridap.RefFEs.NewDivRefFEs._nfaces_array(p,fmoments,cmoments,ft)

# I want to store [Xgp_NF] for all NFS of polytope "Nodes per n-face"
# int_points = Gridap.RefFEs.NewDivRefFEs._nfaces_array(p,fcips,ccips,pt)

nface_moments
nfacedofs = Gridap.RefFEs.NewDivRefFEs._nfacedofs_basis(p,nface_moments)

# Build DOFBasis and RefFE with all this
dof_basis = NewDivDOFBasis(nface_evaluation_points, nface_moments)

divreffe = Gridap.RefFEs.NewDivRefFEs._NewDivRefFE(p,dof_basis,basis,nfacedofs)

##


##
# fns, o = facet_normals(p)
# fns
# o
#
# fns[1]*o[1]

for order in 1:3

  p = Polytope(1,1)

  reffe = RaviartThomasRefFE(p,order)

  dofsb = dofbasis(reffe)

  p = polytope(reffe)

  shb = shfbasis(reffe)

  nfdofs = nfacedofs(reffe)


  b = dofsb
  shfs = reffe.shfbasis
  kk = evaluate(b,shfs)
  @test kk ≈ Matrix(1.0I, size(kk))

  fun(x) = VectorValue(x[1],2*x[1])
  D = 2
  T = VectorValue{2,Float64}
  field = AnalyticalField(fun,D)

  isa(dofsb,DOFBasis{2,T})
  isa(field,Field{2,T})

  evaluate(b,field)
end
##
end # module
