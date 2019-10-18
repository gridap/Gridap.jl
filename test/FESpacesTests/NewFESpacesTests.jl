module NewFEFESpacesTests
##
using Gridap.DOFBases: LagrangianDOFBasis

using Gridap, Test

using Gridap.Helpers

using Base.Cartesian

import Gridap: ∇
# 1D reffe

F = Gridap.NewConformingFESpaces

# Construct the discrete model
n = 2
model = CartesianDiscreteModel(partition=(n,))

diri = [1]
p = Polytope(1,)
order = 2

trian = Triangulation(model)
reffe = LagrangianRefFE(Float64,p,order)

graph = GridGraph(model)
labels = FaceLabels(model)
mydiri_tags = Gridap.BoundaryGrids._setup_tags(labels,diri)

fespace = NewConformingFESpace(reffe,trian,graph,labels,mydiri_tags)
##
# fespace = ConformingFESpace(reffe,trian,graph,labels,[4])
# fespace = MultiShapeConformingFESpace(reffes,trian,graph,labels,mydiri_tags)
# typeof(reffes)
# isa(reffes,IndexCellValue{LagrangianRefFE{1,Float64},1})

# Define manufactured functions
ufun(x) = x[1]*(x[1]-2.0)
ufun_grad(x) = VectorValue(2*x[1]-2.0,)
bfun(x) = -2.0
# ufun(x) = x[1]
# ufun_grad(x) = VectorValue(1.0,)
# bfun(x) = 0.0
∇(::typeof(ufun)) = ufun_grad

# fun = ufun
#
# T = Float64
#
# reffe = fesp.reffe
# dofb = broadcast(dofbasis,reffe)
# trian = fesp.triangulation
# phi = CellGeomap(trian)
# uphys = fun ∘ phi
# celldofs = fesp.cell_eqclass
# nfdofs = fesp.dim_to_nface_eqclass
# E = dof_type(T)
# free_dofs = zeros(E, num_free_dofs(fesp))
# diri_dofs = zeros(E, num_diri_dofs(fesp))
#
# max(collect(broadcast(length,dofb))...)
#
# aux = zeros(E, length(dofb))
# # F._interpolate_values_kernel!(
# # free_dofs,diri_dofs,uphys,celldofs,dofb,aux)
# # function _interpolate_values_kernel!(
# # free_dofs,diri_dofs,uphys,celldofs,dofbs,aux)
#
# dofb[1]
# celldofs[1]
# uphys[1]
#
# evaluate!(dofb[1],uphys[1],aux)
#
#
# b = dofb[1]
# f = uphys[1]
# dofs = aux
#
#
# # function evaluate!(
# # b::LagrangianDOFBasis{D,T},f::Field{D,T},dofs::AbstractVector{E}) where {D,T,E}
# cache = b._cache_field
# evaluate!(f,b.nodes,cache)
# for (node,v) in enumerate(cache)
#   for (comp,vi) in enumerate(v)
#     dof = b.node_and_comp_to_dof[node,comp]
#     dofs[dof] = vi
#   end
# end
# # end
#
# for (dof,imap,l2g) in zip(dofb,uphys,celldofs)
#   evaluate!(dof,imap,aux)
#   # @santiagobadia : Here we should add a method for RT that multiplies by -1
#   # if the face has this cell as the second one in the grid graph
#   # local_to_global_dofs, global_to_local_dofs CellArray
#   for (i,gdof) in enumerate(l2g)
#     if (gdof > 0)
#       free_dofs[gdof] = aux[i]
#     else
#       diri_dofs[-gdof] = aux[i]
#     end
#   end
# end
#
#
#

V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)
##
# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,degree=order*2)

# Define forms
a(v,u) = inner(∇(v), ∇(u))
b(v) = inner(v,bfun)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(a,b,V,U,assem,trian,quad)

# Define the FESolver
ls = BackslashSolver()
solver = LinearFESolver(ls)

# Solve!
uh = solve(solver,op)


# Define exact solution and error
u = CellField(trian,ufun)
e = u - uh

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = a(u,u) + l2(u)

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8

# writevtk(trian,"onedres",cellfields=["uh"=>uh,"e"=>e])
##
end # module
