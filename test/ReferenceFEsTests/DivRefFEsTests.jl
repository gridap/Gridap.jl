module RaviartThomasRefFEsTest

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

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 3

prebasis = QCurlGradMonomialBasis{D}(et,order)

nf_nodes, nf_moments = _RT_nodes_and_moments(et,p,order)

face_own_dofs = _face_own_dofs_from_moments(nf_moments)

dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

v = VectorValue(3.0,0.0)
field = MockField{D}(v)

cache = dof_cache(dof_basis,field)
r = evaluate_dof!(cache, dof_basis, field)

cache = dof_cache(dof_basis,prebasis)
r = evaluate_dof!(cache, dof_basis, prebasis)

display(r)


end # module
