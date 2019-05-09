##
using Numa
using Test

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
##
# @testset ConformingFESpace
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
@test FESpaces.num_free_dofs(fesp) == 9
@test FESpaces.num_fixed_dofs(fesp) == 0
##
fespwd = ConformingFESpace(reffe,trian,gridgr,labels,(1,2,3,4))
@test FESpaces.num_free_dofs(fespwd) == 5
@test FESpaces.num_fixed_dofs(fespwd) == 4

fespwd = ConformingFESpace(reffe,trian,gridgr,labels,(10,))
@test FESpaces.num_free_dofs(fespwd) == 1
@test FESpaces.num_fixed_dofs(fespwd) == 8

@test FESpaces.reffes(fesp) == reffe
@test FESpaces.triangulation(fesp) == trian
@test FESpaces.num_fixed_dofs(fesp) == fesp.num_fixed_dofs
@test FESpaces.num_free_dofs(fesp) == fesp.num_free_dofs
@test FESpaces.nf_eqclass(fesp) == fesp.nf_eqclass
@test FESpaces.cell_eqclass(fesp) == fesp.cell_eqclass
##

fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
funh = interpolate(fun, fesp)

dird = funh.coeffs.gid_to_val_neg.v
funh.coeffs.gid_to_val_pos
##

using Numa.FESpaces: FESpaceWithDirichletData
fespwd = FESpaceWithDirichletData(fesp, dird)
dird
FESpaces.reffes(fespwd)

@test FESpaces.reffes(fespwd) == reffe
@test FESpaces.triangulation(fespwd) == trian
@test FESpaces.num_fixed_dofs(fespwd) == nfixed
@test FESpaces.num_free_dofs(fespwd) == nfree
@test FESpaces.nf_eqclass(fespwd) == dofs
@test FESpaces.cell_eqclass(fespwd) == fesp.cell_eqclass


##
fesphom = FESpaces.TestFESpace(fesp)
@test sum(fesphom.dir_data .== 0.0) == length(fesphom.dir_data)
##

labels
























##################3


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
