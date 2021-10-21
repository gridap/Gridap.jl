using Gridap
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.CellData
using FillArrays

L = 2 # Domain length in each space dimension
D = 2 # Number of spatial dimensions
n = 4 # Partition (i.e., number of cells per space dimension)

pmin = Point(Fill(0,D))
pmax = Point(Fill(L,D))
partition = Tuple(Fill(n,D))
model = simplexify(CartesianDiscreteModel(pmin,pmax,partition))

T = Float64
order = 1
pol = Polytope(Fill(HEX_AXIS,D)...)
reffe = LagrangianRefFE(T,pol,order)

Vₕ = FESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
u(x) = x[1]
Uₕ = TrialFESpace(Vₕ,u)

dv = get_fe_basis(Vₕ)
du = get_trial_fe_basis(Uₕ)

grad_dv = ∇(dv)
grad_du = ∇(du)

grad_dv_array = get_data(grad_dv)
grad_du_array = get_data(grad_du)

Iₖ = lazy_map(Broadcasting(Operation(⋅)),grad_du_array,grad_dv_array)

println("******************")
array_cache(Iₖ)
print_op_tree(Iₖ,showid=true)