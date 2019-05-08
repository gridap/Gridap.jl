##
using Numa, Test

using Numa.Quadratures
using Numa.Polytopes
using Numa.RefFEs
# using Numa.Meshes
using Numa.FESpaces
using Numa.FESpaces: ConformingFESpace

using Numa.Polytopes: PointInt
using Numa.CellValues

using Numa.Maps
using Numa.FieldValues

using Numa.Meshes

using Numa.CellValues: CellVectorByComposition

using Numa.CellMaps

using Numa.Geometry
using Numa.Geometry.Cartesian

import Numa: gradient, ∇

using Numa.CellIntegration
using UnstructuredGrids
##

D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
phi = geomap(trian)

order=1
orders=order*ones(Int64,D)
pol_array = celltypes(trian)
extrusion = PointInt{D}(pol_array[1])
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = ConstantCellValue(basis, ncells(trian))
quad = quadrature(trian,order=2)
gps = Quadratures.coordinates(quad)

fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
# fun(x::Point{2}) = x[1]^2
# gradfun(x::Point{2}) = VectorValue(2*x[1], 0.0)
# gradient(::typeof(fun)) = gradfun

vefcells = veftocells(graph)
is_fixed_vef = zeros(Bool, length(vefcells))
# fesp = ConformingFESpace(reffe,trian,graph,is_fixed_vef)
# assembler = ConformingAssembler(fesp)
# funh = interpolate(fun, fesp)
#
# p_arr = []
# for pi in points(grid)
#   global p_arr
#   p_arr = [p_arr..., pi]
# end
# p_arr
# fun.(p_arr)
# @test funh.coeffs.gid_to_val_pos == fun.(p_arr)
##

# Next, put negative values, etc...
# D=2
# nparts1d = 2
# nparts = nparts1d*ones(Int64,D)
# nparts_t = tuple(nparts...)
# grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
# trian = triangulation(grid) # Generates the Triangulation associated with this grid
# quad = quadrature(trian,order=2)
# fun(x::Point{2}) = 1.0
# ksca = integrate(fun, trian, quad)
# int = sum(ksca)
# @test int ≈ 1.0
##
# @santiagobadia : I would say that this is not the result we want...
# we want the whole integral over the whole domain ... 1.0

# physbasis = attachgeomap(cellb,phi);
# ab(v,u) = inner(∇(v),∇(u)) #+ inner(v,u)
# am(v,u) = inner(v,u)
# V = physbasis;
# U = physbasis;
#
# uphys = fun ∘ phi
#
# evaluate(phi, gps)
# evaluate(uphys, gps)
# ksca = integrate(ab(uphys,uphys),trian,quad)
# typeof(ksca)
# int = sum(ksca)
# @test int ≈ 1.0
##

# kvec = integrate(ab(V,uphys),trian,quad)
# typeof(kvec)
#
# # @santiagobadia : Problem when showing this result, in iterate
# kvec = integrate(am(V,uphys),trian,quad)
# kmat = integrate(ab(V,U),trian,quad)
# # @santiagobadia : Problem in evaluation...
#
# cellvefs = celltovefs(graph)
# vefcells = veftocells(graph)
#
# is_fixed_vef = zeros(Bool, length(vefcells))
# is_fixed_vef[1:4] .= true
#
# fesp = ConformingFESpace(reffe,trian,graph,is_fixed_vef)
# assembler = ConformingAssembler(fesp)
# # assembler.assembly_op_cols
# funh = interpolate(fun, fesp)
# funh.coeffs.gid_to_val_pos
# funh.coeffs.gid_to_val_neg
##

# using Numa.FESpaces: assemble
# sys_vec = assemble(assembler,kvec)
# sys_mat = assemble(assembler,kmat)


####################
#########################33

##
D = 2
model = CartesianDiscreteModel(domain=(0.0,1.0,-1.0,2.0),
        partition=(2,2))

grid = Grid(model,2)
trian = triangulation(grid)
gridgr = FullGridGraph(model)
labels = FaceLabels(model)

order=1; orders=order*ones(Int64,D)
polytope = Polytopes.Polytope(1,1)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
fesp = ConformingFESpace(reffe,trian,gridgr,labels)
dofs, nfree, nfixed = FESpaces.globaldofs(reffe, gridgr::FullGridGraph, labels)
dofs
@test nfree == 1
@test nfixed == 8
##

offset = tuple(length.(fesp.dof_eqclass[1:end])...)
fesp.dof_eqclass

dofs[1]
dofs[2]

using Numa.Geometry: connections
cellvefs = [connections(gridgr,D,i) for i in 0:1]
cellvefs[1]
cellvefs[2]



using Numa.CellValues: IndexCellValueByLocalAppendWithOffset
offset = tuple(length.(fesp.dof_eqclass)...)
cellvefs_dim = [connections(gridgr,D,i) for i in 0:1]
cellvefs = IndexCellValueByLocalAppendWithOffset(offset, cellvefs_dim...)
dofs_all = IndexCellValueByGlobalAppend(dofs...)
dof_eqclass = CellVectorByComposition(cellvefs, dofs_all)

##
