##
using Numa
using Numa.Quadratures
using Numa.CellQuadratures
# using Numa.CellIntegration
# using Numa.CellValues
using Numa.CellFunctions
import Numa.CellIntegration: cellcoordinates, cellbasis

using Numa.CellValues: IndexCellArray
import Numa.CellValues: cellsize
using Numa.Polytopes
using Numa.Polytopes: PointInt
using Numa.RefFEs
using Numa.FieldValues

import Numa: gradient
using Numa.CellValues: ConstantCellValue
include("CellIntegrationTestsMocks.jl")

using Numa.Meshes
using Numa.FESpaces: ConformingFESpace

##
##
D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
order=1
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = CellBasisFromSingleInterpolation(basis)
trian = DummyIntegrationMesh2D(partition=nparts_t)
refquad = TensorProductQuadrature(orders=(2,2))
meshcoords = cellcoordinates(trian)
ncells = length(meshcoords)
quad = ConstantCellQuadrature(refquad,ncells)
phi = geomap(trian)
basis = cellbasis(trian)
physbasis = attachgeomap(basis,phi)
ab(v,u) = inner(∇(v),∇(u)) #+ inner(v,u)
V = physbasis
U = physbasis
# fun(x::Point{2}) = x[1]*x[2] + x[1]
# gradfun(x::Point{2}) = VectorValue(x[2] + 1.0, x[1])
fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
uphys = fun ∘ phi
ksca = integrate(ab(uphys,uphys),trian,quad)
sum(ksca)
kvec = integrate(ab(V,uphys),trian,quad)
typeof(kvec)
kmat = integrate(ab(V,U),trian,quad)
##
# Now let us assemble all these values

# @santiagobadia : I need to harmonize meshes...
mesh = StructHexMesh(nparts)
fesp = ConformingFESpace{D,D,ScalarValue}(reffe,mesh)
using Numa.FESpaces: globaldofs
gldofs = globaldofs(fesp)
ndofs = gldofs[end][end]
using Numa.FESpaces: ConformingAssembler
assembler = ConformingAssembler(fesp)
using SparseArrays
using Numa.FESpaces: Assembler

# @santiagobadia : Not tested yet
using Numa.FESpaces: assemble
sys_vec = assemble(assembler,kvec)
sys_mat = assemble(assembler,kmat)
sys_mat*ones(Float64,16)

fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
f = AnalyticalField(fun,2)
f ∘ phi



using Numa.RefFEs: dofs
dofb = dofs(reffe)
reffe.dofbasis
g = f ∘ phi

using Numa.RefFEs: evaluate
evaluate(dofb,g)
# santiagobadia : How do we want to evaluate fields with geomap ?
