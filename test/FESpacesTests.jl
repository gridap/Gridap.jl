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


#############################################33

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


# dofs[1]
# dofs[2]

# cellvefs = [connections(gridgr,D,i) for i in 0:1]



# POSSIBLE WAYS TO CREATE A "MULTIVALUE CELL ARRAY" FOR CELL TO ALL DOFS,
# SINCE NOW DOFS ARE ALSO INDEXED BY DIM...
#
# struct CellVectorByMultiComposition{T,L<:IndexCellArray{Int,1},V<:IndexCellArray{T,1}} <: IndexCellArray{T,1,CachedArray{T,1,Array{T,1}},1}
#   cell_to_x::Vector{L}
#   x_to_vals::Vector{V}
#   cv::CachedVector{T,Vector{T}}
# end
# # @santiagobadia : For some reason, IndexCellVector{Int} not working
#
# function CellVectorByMultiComposition(cell_to_x::AbstractVector{L},
#   x_to_vals::AbstractVector{V}) where {T,L<:IndexCellArray{Int,1},V<:IndexCellArray{T,1}}
#   # L = eltype(cell_to_x)
#   # V = eltype(x_to_vals)
#   @assert length(cell_to_x) == length(x_to_vals)
#   num_vecs = length(cell_to_x)
#   c = 0
#   for i in 1:num_vecs
#     c += celllength(cell_to_x[i])*celllength(x_to_vals[i])
#   end
#   a = Vector{T}(undef,(c,))
#   cv = CachedArray(a)
#   CellVectorByMultiComposition{T,L,V}(cell_to_x, x_to_vals, cv)
# end
#
# @propagate_inbounds function getindex(self::CellVectorByMultiComposition,cell::Int)
#   setsize!(self.cv,_cellsize(this,cell))
#   l = 1
#   for idim in 1:length(self.cell_to_x)
#     cell_to_x = self.cell_to_x[idim][cell]
#     for x in cell_to_x
#       for val in self.x_to_vals[idim][x]
#         for iv in val
#           self.cv[l] = iv
#           l += 1
#         end
#       end
#     end
#   end
#   self.cv
# end
#
# size(self::CellVectorByMultiComposition) = (length(self.cell_to_x[1]),)
#
# IndexStyle(::Type{CellVectorByMultiComposition{T,L,V}}) where {T,L,V} = IndexLinear()
#
# function cellsize(this::CellVectorByMultiComposition)
#   c = 0
#   for idime in 1:length(this.cell_to_x)
#     c += celllength(this.cell_to_x[idime][i])*celllength(this.x_to_vals[idime][i])
#   end
#   return (c,)
# end
#
# function _cellsize(self::CellVectorByMultiComposition, cell::Int)
#   for idim in 1:length(self.cell_to_x)
#     cell_to_x = self.cell_to_x[idim][cell]
#     l = 0
#     for x in cell_to_x
#       l += length(self.x_to_vals[idim][x])
#     end
#   end
#   return (l,)
# end
#
#
#
# struct CellVectorFromMultiDataAndPtrs{T,V,P} <: IndexCellArray{T,1,CachedSubVector{T,V},1}
#   data::V
#   ptrs::P
#   cv::Vector{CachedSubVector{T,V}}
# end
#
# function CellVectorFromMultiDataAndPtrs(data::AbstractArray{T,1},ptrs::AbstractArray{Int,1}) where T
#   V = typeof(data)
#   P = typeof(ptrs)
#   for i in 1:length(data)
#     cv[i] = CachedSubVector(data[i],0,0)
#   end
#   CellVectorFromMultiDataAndPtrs{T,V,P}(data,ptrs,cv)
# end
#
# @propagate_inbounds function getindex(self::CellVectorFromDataAndPtrs,cell::Int,i::Int)
#   pini = self.ptrs[cell][i]
#   pend = self.ptrs[cell+1][i]-1
#   locate!(self.cv[i],pini,pend)
#   self.cv
# end
#
# size(self::CellVectorFromDataAndPtrs) = (length(self.ptrs)-1,)
#
# IndexStyle(::Type{CellVectorFromDataAndPtrs{T,V,P}}) where {T,V,P} = IndexLinear()
#
# function cellsize(self::CellVectorFromDataAndPtrs)
#   l = 0
#   for i in 1:(length(self.ptrs)-1)
#     c = 0
#     for j in 1:length(self.data)
#       c += self.ptrs[i+1][j]-self.ptrs[i][j]
#     end
#     l = max(l,c)
#   end
#   (l,)
# end
