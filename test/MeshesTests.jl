using Numa, Test

##
spdims = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,spdims)
nprocs1d=1
nprocs = nprocs1d*ones(Int64,spdims)
iproc = tuple(zeros(Int64,spdims)...)
# licells = LexIndexSet(nparts)
# linfs = LexIndexSet(2*ones(Int64,spdims))
polytope = Polytope(ones(Int64,spdims))
@time mesh = StructHexMesh(nparts)
@test mesh.cellvefs[end][mesh.polytope.dimnfs[1]][end] == (nparts1d+1)^spdims
# @test mesh.graph[mesh.polytope.dimfs[spdims+1],:][end] == 25#nparts1d^spdims
##
