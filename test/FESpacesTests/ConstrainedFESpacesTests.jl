module ConstrainedFESpacesTests

using Test
using Gridap
using Gridap.ConstrainedFESpaces: VectorOfTwoParts

D = 2
dom = fill(4,D)
model = CartesianDiscreteModel(partition=tuple(dom...))

order = 1
diritag = "boundary"
_fespace = CLagrangianFESpace(Float64,model,order,diritag)

fixeddofs = [2,5]
fespace = ConstrainedFESpace(_fespace,fixeddofs)

ufun(x) = x[1] + x[2] + 2.0

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

bh = FEBasis(fespace)
uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

nfree = 7

ndiri = 16

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

g(x) = 2.0

U = TrialFESpace(fespace,g)

u(x) = x[1] + x[2] + 2.0

uh = interpolate(fespace,u)

#writevtk(trian,"trian",cellfields=["uh"=>uh])

nfree = 10
ndiri = 4
fixeddofs = [3,1,5,7]
nfixed = length(fixeddofs)
nfree_new = nfree - nfixed
dof_to_new_dof = zeros(Int,nfree)
is_fixed = fill(false,nfree)
is_fixed[fixeddofs] .= true
is_free = is_fixed.==false
dof_to_new_dof[is_free] .= 1:nfree_new
dof_to_new_dof[is_fixed] .= (-(ndiri+1)):-1:(-(ndiri+nfixed))

freevals=[20,40,60,80,90,100]
fixedvals=[30,10,50,70]

v = VectorOfTwoParts(dof_to_new_dof,freevals,fixedvals,ndiri)
@test all(v[is_free] .== freevals)
@test all(v[is_fixed] .== fixedvals)

end # module
