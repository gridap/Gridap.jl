# module CDLagrangianRefFEsTests

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs

# T = VectorValue{2,Float64}
T = Float64
p = SEGMENT
orders = (2,)
cont = (DISC,)
D = 1

# lrfe = LagrangianRefFE(T,p,orders)
# grfe = lrfe.data
#
# X = Gridap.ReferenceFEs
# face_own_nodes = X._compute_cd_face_own_nodes(p,orders,cont)
#
# face_own_dofs = X._generate_face_own_dofs(face_own_nodes, grfe.dofs.node_and_comp_to_dof)
#
# reffaces = X.compute_lagrangian_reffaces(T,p,orders)
#
# _reffaces = vcat(reffaces...)
#
# ndofsf = length(grfe.dofs.dof_to_node)
#
# face_own_dofs_permutations = X._trivial_face_own_dofs_permutations(face_own_dofs)
#
# face_dofs = X._generate_face_dofs(ndofs,face_own_dofs,p,_reffaces)
# # face_dofs = Vector{Vector{Int}}()
#
# GenericRefFE(
#                 grfe.ndofs,
#                 grfe.polytope,
#                 grfe.prebasis,
#                 grfe.dofs,
#                 grfe.face_own_dofs,
#                 grfe.face_own_dofs_permutations,
#                 grfe.face_dofs,
#                 grfe.shapefuns)
#

reffe = CDLagrangianRefFE(T,SEGMENT,(2,),(DISC,))
test_reference_fe(reffe)

reffe = CDLagrangianRefFE(T,QUAD,(2,2),(CONT,DISC))
test_reference_fe(reffe)

T = VectorValue{2,Float64}
reffe = CDLagrangianRefFE(T,HEX,(2,2,2),(CONT,CONT,DISC))
test_reference_fe(reffe)

end # module
