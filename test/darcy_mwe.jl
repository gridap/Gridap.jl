using Gridap
using Gridap.Geometry
using Gridap.FESpaces

D = 2
order = 2
domain = ifelse(D==2, (0,1,0,1), (0,1,0,1,0,1))
partition = Tuple(fill(3, D))
model = simplexify(CartesianDiscreteModel(domain,partition); positive=true)

Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*(order+1))

# u_exact: degree-2, div-free, in [P2]³ ⊂ RT2.
#   Derived from stream potential A=(0,0,x²y): curl(A) = (x², -2xy, 0), div=0.
# p_exact: degree-2, in P2.
# Body force g = u_exact + ∇p_exact makes (u_exact,p_exact) the exact Darcy solution.
u_exact(x) = ifelse(D==2,VectorValue(x[1]^2, -2*x[1]*x[2]), VectorValue(x[1]^2, -2*x[1]*x[2], 0.0))
p_tilde(x) = x[1]^2 + x[2]^2
p_mean = sum(∫(p_tilde)dΩ) / sum(∫(1)dΩ)
p_exact(x) = p_tilde(x) - p_mean

f(x) = u_exact(x) - ∇(p_exact)(x)

reffe_u = ReferenceFE(raviart_thomas, Float64, order; change_dof=false)
reffe_p = ReferenceFE(lagrangian, Float64, order)

V = FESpace(model, reffe_u; conformity=:HDiv, dirichlet_tags="boundary")
Q = FESpace(model, reffe_p; conformity=:L2, constraint=:zeromean)
U = TrialFESpace(V, u_exact)
X = MultiFieldFESpace([U, Q])
Y = MultiFieldFESpace([V, Q])

a((u,p),(v,q)) = ∫( u⋅v + p*(∇⋅v) + (∇⋅u)*q )dΩ
l((v,q))       = ∫( f⋅v )dΩ

op = AffineFEOperator(a, l, X, Y)
xh = solve(op)
uh, ph = xh

eu = uh - u_exact
ep = ph - p_exact
println("‖u - u_exact‖_L2 = ", sqrt(sum(∫(eu⋅eu)dΩ)))
println("‖p - p_exact‖_L2 = ", sqrt(sum(∫(ep*ep)dΩ)))
