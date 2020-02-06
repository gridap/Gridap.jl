# get_faces restricted to dim enumeration seems useless...
##

using Gridap
using Gridap.ReferenceFEs
using Gridap.Helpers
using Gridap.Polynomials
using Gridap.Fields
using Gridap.Integration
using Gridap.Fields: MockField, MockBasis, OtherMockBasis
using Gridap.Arrays

# I cannot include it in ReferenceFEs because of a cycle...
include("../../src/ReferenceFEs/MomentReferenceFEs.jl")

p = Polytope(HEX_AXIS, HEX_AXIS)
D = num_dims(p)
et = Float64
order = 3


# Given p, D, et, order, we must create the RT constructor, code below...

if !(all(get_extrusion(p).array .== HEX_AXIS))
  @notimplemented
end
# 1. Prebasis
prebasis = QCurlGradMonomialBasis{D}(et,order)

# Nface nodes, moments, and prebasis evaluated at nodes
nf_nodes, nf_moments, pb_moments = _initialize_arrays(prebasis,p)

# Face values
fcips, fmoments = _RT_face_values(p,et,order)
nf_nodes,nf_moments,pb_moments = _insert_nface_values!(nf_nodes,nf_moments,pb_moments,prebasis,fcips,fmoments,p,dim(p)-1)

# Cell values
if (order > 1)

  ccips, cmoments = _RT_cell_values(p,et,order)
  nf_nodes,nf_moments,pb_moments = _insert_nface_values!(nf_nodes,nf_moments,pb_moments,prebasis,ccips,cmoments,p,dim(p))

end

# Change of basis matrix, inv([DF,DC])
cob = inv(hcat(pb_moments))
basis = change_basis(prebasis,cob)

nfacedofs = _nfacedofs_basis(p,nf_moments)

# Build DOFBasis and RefFE with all this
dof_basis = _GenericDofBasis(nf_nodes, nf_moments)

# This part is missing, it is just to create the struct and implement the API
# divreffe = _GenericRefFE(p,dof_basis,basis,nfacedofs)
# return divreffe

# Here some tests I have been doing
b = dof_basis

v = VectorValue(3.0,0.0)
d = 2
field = MockField{d}(v)

cache = dof_cache(dof_basis,field)
evaluate_dof!(cache, b, field)

ndof = 8
w = fill(v,ndof)
f = OtherMockBasis{d}(ndof)
basis = MockBasis{d}(v,ndof)

cache = dof_cache(dof_basis,basis)
evaluate_dof!(cache, b, basis)

# - [`dof_cache`](@ref)
# - [`evaluate_dof!`](@ref)
