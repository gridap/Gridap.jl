module DarcyFEOperatorsTests

##
using Test
using Gridap
import Gridap: ∇

import LinearAlgebra: det
using Base: transpose
##

# model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(2,2))
# model.cgrid
# ps = points(model.cgrid)
# typeof(ps)
# collect(ps)
#
#
#
# f(x) = Point(x[1],5*x[2])
#
# broadcast(f,ps)
#
#
# struct MappedCartesianGrid <: CartesianGrid
#   map
#
# end
# with mapping
#
#
# trian = Triangulation(model)
#
#
# phi = CellGeomap(trian)
# jac = gradient(phi)
# jact = transpose(jac)
#
#
# piola_map = inv(jact)
#
#
# graph = GridGraph(model)
# pols = CellPolytopes(Grid(model))
# pt = pols.value
#
# order = 1
# _reffe = NedelecRefFE(pt,Float64,order)
#
# _labels = FaceLabels(model)
#
# # V = ConformingFESpace(_reffe,trian,graph,_labels,[5,6])
# V = CurlConformingFESpace(_reffe,trian,graph,_labels,[5,6])
##

u(x) = VectorValue(x[1]*x[2], -0.5*x[2]^2)
∇u(x) = TensorValue(x[2],0.0,x[1],-x[2])
divu(x) = 0.0
∇(::typeof(u)) = ∇u

const kinv_1 = TensorValue(1.0,0.0,0.0,1.0)
r(x) = kinv_1*u(x)

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(1,1))
order = 2
V = FESpace( reffe=:RaviartThomas, conformity=:HDiv, order=order, model=model, diritags = [5,6])

V_0 = TestFESpace(V)
V_g = TrialFESpace(V,u)

trian = Triangulation(model)
quad = CellQuadrature(trian,degree=order*2)

@law σ(x,u) = kinv_1*u

function a(v,u)
  inner(v,σ(u))
end

function b_Ω(v)
  inner(v,r)
end

t_Ω = AffineFETerm(a,b_Ω,trian,quad)

# op = LinearFEOperator(V_0,V_g,t_Ω)

w = weights(quad)
z = coordinates(quad)

evaluate(V.cellbasis,z)


phi = CellGeomap(trian)
jac = gradient(phi)
detjac = det(jac)
#
piola_map = (1.0/detjac)*jac
s = V.cellbasis
pbasis = Gridap.CellFieldsOperations._merge_val_and_grad(piola_map*s.val,piola_map*s.grad)
#
vpm = evaluate(piola_map,z)
size(vpm[1])
vs = evaluate(kinv_1*s,z)
aaa = kinv_1*s
typeof(aaa)
length(vs[1])
@show vs[1]
@show vpm[1]
evaluate(piola_map.*transpose(s.val),z)

(vpm[1].*vs[1]')'
length(z[1])
evaluate(s,z)

uh = solve(op)

e1 = u - uh

# Define norms to measure the error
l2(u) = inner(u,u)

# Compute errors
e1l2 = sqrt(sum( integrate(l2(e1),trian,quad) ))

# writevtk(trian,"/home/santiago/github-repos/Gridap/tmp/darcyresults",
         # cellfields=["uh"=>uh,"ph"=>ph, "euh"=> e1, "eph"=> e2])

@test e1l2 < 1.e-8
##
end # module
