include("../../src/FESpaces/FETerms.jl")
module FETermsTests

using Test
using Gridap
using ..FETerms

using ..FETerms: _restrict_if_needed

ufun(x) = x[1] + x[2]
bfun(x) = x[1]
gfun(x) = x[2]

model = CartesianDiscreteModel(partition=(4,4))

order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

tags = [7,6]
btrian = BoundaryTriangulation(model,tags)
bquad = CellQuadrature(btrian,order=2)

v = FEBasis(V)
du = FEBasis(U)
uhd = zero(U)
uh = interpolate(U,ufun)

bfield = CellField(trian,bfun)
gfield = CellField(btrian,gfun)

a(v,u) = inner(v,u)
l(v) = inner(v,bfield)
g(v) = inner(v,gfield)

t1 = AffineFETerm(a,l,trian,quad)

cm = setup_cell_matrix(t1,v,du)
@test isa(cm,CellMatrix)
@test length(cm) == ncells(trian)

cv = setup_cell_vector(t1,v,uhd)
@test isa(cv,CellVector)
@test length(cv) == ncells(trian)

cn = setup_cell_ids(t1)
@test isa(cn,CellNumber)
@test length(cn) == ncells(trian)

t2 = AffineFETerm(a,g,btrian,bquad)

cm = setup_cell_matrix(t2,v,du)

#_v = _restrict_if_needed(v,btrian)
#_u = _restrict_if_needed(du,btrian)
#
#@show isa(_v,FEBasis)
#@show isa(_u,FEBasis)
#cm =  integrate(a(_v,_u),btrian,bquad)
#
#q = coordinates(bquad)
#w = weights(bquad)
#
#cmap = a(_v,_u)
#cmq = evaluate(cmap,q)
#
##@show cmq
#
#phi = CellGeomap(btrian)
#j = evaluate(∇(phi),q)
#@show size(cmq[1])
#
#@show size(meas(j)[1])
#
#cmq[1] .* meas(j)[1]

end #module
