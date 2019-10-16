module MappedCartesianGeometryTests

using Test
using Gridap


###
F = Gridap.CartesianGeometry

model1 = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0),partition=(1,1))
points(model1.cgrid)
model1.cgrid
# trian1 = Triangulation(model1)
# points(trian1.grid)
# phi1 = CellGeomap(trian1)

gmap(x) = 2*x
model2 = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0),partition=(1,1),map=gmap)
points(model2.cgrid)
trian = Triangulation(model2)
# points(trian2.grid)
phi = CellGeomap(trian)

##

# CellPoints(trian1)

# CellPoints(trian2)


##




quad = CellQuadrature(trian,degree=2)
z = coordinates(quad)
##
self = trian1
Gridap.CellValuesGallery.CellVectorFromLocalToGlobal(cells(self.grid),points(self.grid))

self = trian2
cells(trian1.grid)

##
# phi1.val.cellvalues
# typeof(phi1.val.cellvalues[2]) == typeof(phi2.val.cellvalues[2])
#
# phi1.val.cellvalues[2].lid_to_gid
# phi1.val.cellvalues[2].gid_to_val
# phi1.val.cellvalues[2].cv

  # gid_to_val::V
  # cv::CachedArray{T,1,Array{T,1}}
# phi2.val.cellvalues[2].lid_to_gid
# phi2.val.cellvalues[2].gid_to_val[1]
# phi2.val.cellvalues[2].cv
##














# getindex(phi2.val,1)
# vals = Gridap.CellMapApply._getvalues(1,phi2.val.cellvalues...)
# self = phi1
# vals = Gridap.CellMapApply._getvalues(1,self.cellvalues...)
#
#
#
evaluate(phi1[1],z[1])
evaluate(phi2[1],z[1])





# evaluateI
# pval = collect(phi.val)
# pval[1].inputs
#
# isa(phi,Gridap.CellFieldsOperations.IndexCellFieldLikeAndGradient)
# typeof(phi)
#
# ##
# phi
# uphys = fun ∘ phi
#
# celldofs = fesp.cell_eqclass
# nfdofs = fesp.dim_to_nface_eqclass
# E = dof_type(T)
# free_dofs = zeros(E, num_free_dofs(fesp))
# diri_dofs = zeros(E, num_diri_dofs(fesp))
# aux = zeros(E, length(dofb))
# _interpolate_values_kernel!(
#
# # F._compute_points(cg)
#
# cg[1]
#
# # Here we have the points
# model.ugrid.points

# model.gridgraph.data

end # module
