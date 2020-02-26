module StateLawsTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Integration
using Gridap.FESpaces
using Gridap.Fields

domain = (0,1,0,1)
n = 3
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

order = 3
V = TestFESpace(model=model,reffe=:Lagrangian,valuetype=Float64,order=order,conformity=:L2)
U = TrialFESpace(V)

v = get_cell_basis(V)
du = get_cell_basis(V)
uh = FEFunction(U,rand(num_free_dofs(U)))

@statelaw function foo(u,s,r)
  w = u
  w, s, r+1
end

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

s = QPointCellField(0.0,trian,quad)
r = QPointCellField(0.0,trian,quad)

wh = foo(uh,s,r)

q = get_coordinates(quad)
wh_q = evaluate(wh,q)
uh_q = evaluate(uh,q)


function loop(a)
  cache = array_cache(a)
  for k in 1:length(a)
    ak = getindex!(cache,a,k)
  end
end

r_q = evaluate(r,q)
@test r_q[end] == zeros(size(q[end]))

loop(wh_q)

@test r_q[end] == ones(size(q[end]))

end # module
