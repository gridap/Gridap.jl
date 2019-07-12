include("../../src/Integration/SkeletonCellFields.jl")
module SkeletonCellFieldsTests

using Test
using Gridap
import Gridap: ∇
using ..SkeletonCellFields

model = CartesianDiscreteModel(partition=(3,4,2))

order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

trian = Triangulation(model)

tags = [27,]
strian = SkeletonTriangulation(model,tags)
squad = CellQuadrature(strian,order=2)
s = coordinates(squad)

ufun(x) = 2.0*x[1]^2 * x[2]^2
grad_ufun(x) = VectorValue(4.0*x[1]*x[2]^2, 4.0*x[1]^2*x[2])
∇(::typeof(ufun)) = grad_ufun

u = CellField(trian,ufun)

uh = interpolate(fespace,ufun)

su = restrict(u,strian)

suh = restrict(uh,strian)

jsu = jump(su)
jsu_s = evaluate(jsu,s)

jsuh = jump(suh)
jsuh_s = evaluate(jsuh,s)

jsu_grad = jump(∇(su))
jsu_grad_s = evaluate(jsu_grad,s)

jsuh_grad = jump(∇(suh))
jsuh_grad_s = evaluate(jsuh_grad,s)

writevtk(strian,"strian",cellfields=[
 "jumpu"=>jsu,"jumpu_grad"=>jsu_grad,
 "jumph"=>jsuh,"jumpuh_grad"=>jsuh_grad ])

end # module
