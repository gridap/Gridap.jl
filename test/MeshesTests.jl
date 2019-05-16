using Gridap, Test

using Gridap.Polytopes
using Gridap.Polytopes: PointInt
using Gridap.Meshes
using Base.Cartesian

##
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nprocs1d=1
nprocs = nprocs1d*ones(Int64,D)
iproc = tuple(zeros(Int64,D)...)
# licells = LexIndexSet(nparts)
# linfs = LexIndexSet(2*ones(Int64,D))
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytope(extrusion)
@time mesh = Meshes.StructHexMesh(nparts)
# @test mesh.cellvefs[end][mesh.polytope.dimnfs[1]][end] == (nparts1d+1)^D
# @test mesh.graph[mesh.polytope.dimfs[D+1],:][end] == 25#nparts1d^D
##
