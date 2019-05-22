##
using Gridap
using Test

using Gridap.Quadratures
using Gridap.Polytopes
using Gridap.Polytopes: PointInt
using Gridap.RefFEs
using Gridap.FESpaces
using Gridap.FESpaces: ConformingFESpace

using Gridap.FieldValues
using Gridap.Maps
using Gridap.CellValues
using Gridap.CellValues.ConstantCellValues
using Gridap.CellValues.Wrappers
using Gridap.CellMaps

using Gridap.Geometry
using Gridap.Geometry.Cartesian

import Gridap: gradient, ∇

using Gridap.CellIntegration
using UnstructuredGrids
##
abstract type mystyle end
struct style1 <: mystyle end
struct style2 <: mystyle end

abstract type str{T} end
mystyle(::str)::mystyle = error("Not defined")

struct str1{T} <: str{T}
  a::T
end
mystyle(::Type{str1{T}}) where {T}= style1

struct str2{T} <: str{T}
  a::T
end
mystyle(::Type{str2{T}}) where {T} = style2

struct strc{T,S<:str{T}}
  a::T
  s::S
end
mystyle(::Type{strc{T,S}}) where {T,S<:str1{T}} = style1
mystyle(::Type{strc{T,S}}) where {T,S<:str2{T}} = style2

s1 = str1(10)
s2 = str2(10)
s1c = strc(20,s1)
s2c = strc(20,s2)



mystyle(typeof(s1))
mystyle(typeof(s2))
mystyle(typeof(s1c))
mystyle(typeof(s2c))




##
D = 2
model = CartesianDiscreteModel(domain=(0.0,1.0,-1.0,2.0),
        partition=(2,2))

grid = Grid(model,2)
trian = Triangulation(grid)
gridgr = FullGridGraph(model)
labels = FaceLabels(model)
dtags = (1,2,3,4)

order=1; orders=order*ones(Int64,D)
polytope = Polytopes.Polytope(1,1)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
##
# @testset ConformingFESpace
fesp = ConformingFESpace(reffe,trian,gridgr,labels,dtags)
D = 2
model = CartesianDiscreteModel(domain=(0.0,1.0,-1.0,2.0),
        partition=(2,2))

grid = Grid(model,2)
trian = Triangulation(grid)
gridgr = FullGridGraph(model)
labels = FaceLabels(model)

order=1; orders=order*ones(Int64,D)
polytope = Polytopes.Polytope(1,1)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
fesp = ConformingFESpace(reffe,trian,gridgr,labels)
@test FESpaces.num_free_dofs(fesp) == 9
@test FESpaces.num_fixed_dofs(fesp) == 0
@test FESpaces.MeshConformity(fesp) == FESpaces.ConformingMesh()

##
fesp = ConformingFESpace(reffe,trian,gridgr,labels,(10,))
@test FESpaces.num_free_dofs(fesp) == 1
@test FESpaces.num_fixed_dofs(fesp) == 8

fesp = ConformingFESpace(reffe,trian,gridgr,labels,(1,2,3,4))
@test FESpaces.num_free_dofs(fesp) == 5
@test FESpaces.num_fixed_dofs(fesp) == 4

@test FESpaces._reffes(fesp) == reffe
@test FESpaces._triangulation(fesp) == trian
@test FESpaces._gridgraph(fesp) == fesp._gridgraph
@test FESpaces._labels(fesp) == fesp._labels
@test FESpaces.num_fixed_dofs(fesp) == fesp.num_fixed_dofs
@test FESpaces.num_free_dofs(fesp) == fesp.num_free_dofs
@test FESpaces.nf_dofs(fesp) == fesp.nf_dofs
@test FESpaces.cell_eqclass(fesp) == fesp.cell_eqclass
##
# @testset FESpaceWithDirichletData
using Gridap.FESpaces: FESpaceWithDirichletData
dird = zeros(Float64, FESpaces.num_fixed_dofs(fesp))
fespwd = FESpaceWithDirichletData(fesp, dird)

@test FESpaces._reffes(fespwd) == reffe
@test FESpaces._gridgraph(fespwd) == fesp._gridgraph
@test FESpaces._labels(fespwd) == fesp._labels
@test FESpaces._triangulation(fespwd) == trian
@test FESpaces.num_fixed_dofs(fespwd) == fesp.num_fixed_dofs
@test FESpaces.num_free_dofs(fespwd) == fesp.num_free_dofs
@test FESpaces.nf_dofs(fespwd) == fesp.nf_dofs
@test FESpaces.cell_eqclass(fespwd) == fesp.cell_eqclass

fesphom = FESpaces.TestFESpace(fesp)
@test sum(fesphom.dir_data .== 0.0) == length(fesphom.dir_data)
# end
##

# TrialFESpace
##
fun1(x::Point{2}) = 5.0
fun2(x::Point{2}) = 4.0
fun3(x::Point{2}) = 3.0
func = [fun1, fun2, fun3, fun2]
fesp = ConformingFESpace(reffe,trian,gridgr,labels,dtags)

fesp.nf_dofs
nf_dofs_all = FESpaces.nf_dofs(fesp)
dtags = FESpaces.dir_tags(fesp)

FESpaces._is_fixed(2,(),labels)

fixed_dofs = FESpaces.interpolate_dirichlet_data(func, fesp)
FESpaces.MeshConformity(fesp)
##
# FEFunction and Interpolate
fesphom = FESpaces.TestFESpace(fesp)
@test FESpaces.MeshConformity(typeof(fesphom)) == FESpaces.ConformingMesh()
fespwnhdd= FESpaces.TrialFESpace(fesp, func, labels)
@test FESpaces.MeshConformity(typeof(fespwnhdd)) == FESpaces.ConformingMesh()
@test fespwnhdd.dir_data == fixed_dofs
fh0 = FESpaces.interpolate(fun1, fesp)
fh1 = FESpaces.interpolate(fun1, fespwnhdd)
fh2 = FESpaces.interpolate(fun1, fesphom)
sum(FESpaces.fixed_dofs(fh2).== 0) == 4
@test FESpaces.fixed_dofs(fh0) != FESpaces.fixed_dofs(fh1)
@test FESpaces.free_dofs(fh0) == FESpaces.free_dofs(fh1) == FESpaces.free_dofs(fh2)
##






















##################3


D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = Triangulation(grid) # Generates the Triangulation associated with this grid
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
quad = CellQuadrature(trian,order=2)
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
# trian = Triangulation(grid) # Generates the Triangulation associated with this grid
# quad = CellQuadrature(trian,order=2)
# fun(x::Point{2}) = 1.0
# ksca = integrate(fun, trian, quad)
# int = sum(ksca)
# @test int ≈ 1.0
##
# @santiagobadia : I would say that this is not the result we want...
# we want the whole integral over the whole domain ... 1.0

# physbasis = attachgeomap(cellb,phi);
# ab(v,u) = varinner(∇(v),∇(u)) #+ varinner(v,u)
# am(v,u) = varinner(v,u)
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

# using Gridap.FESpaces: assemble
# sys_vec = assemble(assembler,kvec)
# sys_mat = assemble(assembler,kmat)


####################
#########################33
