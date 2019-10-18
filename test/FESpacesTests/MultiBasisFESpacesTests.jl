module MultiBasisFESpaces
##
using Gridap.DOFBases: LagrangianDOFBasis

using Gridap, Test

using Gridap.Helpers

using Base.Cartesian

import Gridap: ∇
# 1D reffe

include("ScaledNodesLagrangianRefFEs.jl")
include("MultiBasisFESpaces.jl")
F = Gridap.ConformingFESpaces

D = 2
n = 1
# ℓ = 1.0
β = 400.0
# h = 1.0 / n
# h0 = 1.0 / (n-1+β)
# r = h/h0
# e = h/(ℓ-h0*(n-1))
# ℓ0 = n*h
# geomap = ( x -> ( x[1] <= (n-1)*h ) ? x[1]/r : (n-1)*h0+(x[1]-(n-1)*h)/e )
# geomap = ( x -> x[1] <= (ℓ-h) ? x[1] : (x[1]-ℓ+h)*β+ℓ-h)
# newl = geomap(1.0)
# typeof(newl)
geomap = ( x -> x )
model = CartesianDiscreteModel(partition=(n,n),map=geomap)

# Construct the discrete model
ufun(x) = x[1]*(x[1]-2.0)
ufun_grad(x) = VectorValue(2*x[1]-2.0,0.0)
bfun(x) = -2.0
∇(::typeof(ufun)) = ufun_grad

function solve_problem(β,model,ufun)

  # diri = [7]
  diri = [1,2,3,4,5,6,7,8]
  diri = "boundary"
  p = Polytope(1,1)
  order = 2

  trian = Triangulation(model)
  reffe = LagrangianRefFE(Float64,p,order)
  sreffe = ScaledNodesLagrangianRefFE(Float64,p,order,1.0/β)
  values = LagrangianRefFE{2,Float64}[reffe,sreffe]
  pointer = [ i == ncells(trian) ? 2 : 1 for i in 1:ncells(trian)]
  reffes = CompressedCellValue(values,pointer)
  typeof(reffes)

  graph = GridGraph(model)
  labels = FaceLabels(model)
  mydiri_tags = Gridap.BoundaryGrids._setup_tags(labels,diri)

  fespace = MultiShapeConformingFESpace(reffes,trian,graph,labels,mydiri_tags)

  V = TestFESpace(fespace)
  U = TrialFESpace(fespace,ufun)
  trian = Triangulation(model)
  quad = CellQuadrature(trian,degree=order*2)

  a(v,u) = inner(∇(v), ∇(u))
  b(v) = inner(v,bfun)

  assem = SparseMatrixAssembler(V,U)
  op = LinearFEOperator(a,b,V,U,assem,trian,quad)
  ls = BackslashSolver()
  solver = LinearFESolver(ls)
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

  @show el2 #< 1.e-8
  @show eh1 #< 1.e-8

  return op
end


op = solve_problem(1,model,ufun)
op.mat

op2 = solve_problem(100,model,ufun)
op2.mat

# writevtk(trian,"onedres",cellfields=["uh"=>uh,"e"=>e])
##
end # module
