module SkeletonCellFieldsTests

using Test
using Gridap
import Gridap: ∇

model = CartesianDiscreteModel(partition=(3,4,2))

order = 1
diritag = "boundary"
fespace = H1ConformingFESpace(Float64,model,order,diritag)

trian = Triangulation(model)

tags = [27,]
strian = SkeletonTriangulation(model,tags)
squad = CellQuadrature(strian,degree=2)
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

msu = mean(su)
msu_s = evaluate(msu,s)

jsuh = jump(suh)
jsuh_s = evaluate(jsuh,s)

msuh = mean(suh)
msuh_s = evaluate(msuh,s)

jsu_grad = jump(∇(su))
jsu_grad_s = evaluate(jsu_grad,s)

msu_grad = mean(∇(su))
msu_grad_s = evaluate(msu_grad,s)

jsuh_grad = jump(∇(suh))
jsuh_grad_s = evaluate(jsuh_grad,s)

msuh_grad = mean(∇(suh))
msuh_grad_s = evaluate(msuh_grad,s)

#writevtk(strian,"strian",cellfields=[
# "jumpu"=>jsu,"meanu"=>msu,"jumpu_grad"=>jsu_grad,"meanu_grad"=>msu_grad,
# "jumph"=>jsuh,"meanh"=>msuh,"jumpuh_grad"=>jsuh_grad,"meanh_grad"=>msuh_grad ])

v = CellBasis(trian)
du = v

u_Γ = restrict(u,strian)
du_Γ = restrict(v,strian)
v_Γ = restrict(v,strian)

cm = inner(jump(v_Γ),jump(u_Γ))

cv = integrate(cm,strian,squad)
_ = collect(cv.cellvector1)
_ = collect(cv.cellvector2)

cm = inner(mean(v_Γ),jump(u_Γ))

cv = integrate(cm,strian,squad)
_ = collect(cv.cellvector1)
_ = collect(cv.cellvector2)

cm = inner(jump(v_Γ),jump(du_Γ))

cv = integrate(cm,strian,squad)
_ = collect(cv.cellmatrix11)
_ = collect(cv.cellmatrix12)
_ = collect(cv.cellmatrix21)
_ = collect(cv.cellmatrix22)

cm = inner(mean(v_Γ),jump(du_Γ))

cv = integrate(cm,strian,squad)
_ = collect(cv.cellmatrix11)
_ = collect(cv.cellmatrix12)
_ = collect(cv.cellmatrix21)
_ = collect(cv.cellmatrix22)

cm = inner(jump(u_Γ),jump(u_Γ))
cv = integrate(cm,strian,squad)
_ = collect(cv)

phi = CellGeomap(trian)
phi_Γ = restrict(phi,strian)

x = evaluate(jump(phi_Γ),s)
z = Point(0.0,0.0,0.0)
zs = fill(z,4)
@test all([ xi .+ 1 ≈ zs .+ 1 for xi in x ])

x = evaluate(mean(phi_Γ),s)
x2 = evaluate(phi_Γ.cellfield1,s)
@test x ≈ x2

end # module
