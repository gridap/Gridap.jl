module BoundaryCellFieldsTests

using Test
using Gridap
using Gridap.CellArrayApply: IndexCellArrayFromKernel
import Gridap: ∇

tags = [23,24,25]
model = CartesianDiscreteModel(partition=(3,4,2))

trian = Triangulation(model)
btrian = BoundaryTriangulation(model,tags)

ufun(x) = 2*x[1] + x[2]
grad_ufun(x) = VectorValue(2.0,1.0)
∇(::typeof(ufun)) = grad_ufun

u = CellField(trian,ufun)

bu = restrict(u,btrian)
bu_grad = ∇(bu)

bquad = CellQuadrature(btrian,order=0)
q = coordinates(bquad)

quad = CellQuadrature(trian,order=2)
s = coordinates(quad)

bphi = CellGeomap(btrian)

x = evaluate(bphi,q)

buq = evaluate(bu,q)
@test isa(buq,IndexCellArrayFromKernel)
bu_grad_q = evaluate(bu_grad,q)

phi = CellGeomap(trian)
bphi2 = restrict(phi,btrian)

x2 = evaluate(bphi2,q)

@test x ≈ x2

#writevtk(btrian,"btrian",cellfields=["u"=>bu,"u_grad"=>bu_grad])

#writevtk(x,"x",pointdata=["u"=>buq])

cn = integrate(bu,btrian,bquad)
@test isa(cn,CellNumber)
_ = collect(cn)

v = CellBasis(trian)
bv = restrict(v,btrian)
cm = varinner(bv,bu)

ca = integrate(cm,btrian,bquad)
@test isa(ca,CellArray{Float64,1})
_ = collect(ca)

cm = varinner(bv,bv)
ca = integrate(cm,btrian,bquad)
@test isa(ca,CellArray{Float64,2})
_ = collect(ca)

end # module
