module CrouzeixRaviartTests

using Test
using Gridap
using Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces, Gridap.Arrays, Gridap.TensorValues
using Gridap.Helpers

# Only defined for simplices and order 1
@test_throws ErrorException CrouzeixRaviartRefFE(Float64,QUAD,1)
@test_throws ErrorException CrouzeixRaviartRefFE(Float64,TRI,2)

reffe = ReferenceFE(TRI, crouzeix_raviart, 1)
reffec = CrouzeixRaviartRefFE(Float64, TRI, 1)
@test reffe == reffec
test_reference_fe(reffe)
@test Conformity(reffe, :L2) == L2Conformity()
@test_throws ErrorException Conformity(reffe, :H1)

reffe = ReferenceFE(TET, crouzeix_raviart, 1)
reffec = CrouzeixRaviartRefFE(Float64, TET, 1)
@test reffe == reffec
test_reference_fe(reffe)
@test Conformity(reffe, :L2) == L2Conformity()
@test_throws ErrorException Conformity(reffe, :H1)

function solve_crScalarPoisson(partition, cells, u_exact)
  f(x) = - Δ(u_exact)(x)

  model = simplexify(CartesianDiscreteModel(partition, cells))

  # reffe = CrouzeixRaviartRefFE(VectorValue{2,Float64},TRI,1)
  reffe = CrouzeixRaviartRefFE(Float64,TRI,1)
  V = FESpace(model,reffe,dirichlet_tags="boundary")
  U = TrialFESpace(V,u_exact)

  Ω = Triangulation(model)
  dΩ = Measure(Ω,3)

  a(u,v) = ∫( ∇(u)⋅∇(v) )*dΩ
  l(v) = ∫( f*v )*dΩ

  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)
  e = uh-u_exact

  return sqrt(sum( ∫( e⋅e )*dΩ ))
end


function solve_crVectorPoisson(partition, cells, u_exact)
  f(x) = - Δ(u_exact)(x)

  model = simplexify(CartesianDiscreteModel(partition, cells))

  reffe = CrouzeixRaviartRefFE(VectorValue{2,Float64},TRI,1)
  V = FESpace(model,reffe,dirichlet_tags="boundary")
  U = TrialFESpace(V, u_exact)

  Ω = Triangulation(model)
  dΩ = Measure(Ω,3)

  a(u,v) = ∫( ∇(u) ⊙ ∇(v) )*dΩ
  l(v) = ∫( f ⋅ v )*dΩ

  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)
  e = uh-u_exact

  return sqrt(sum( ∫( e⋅e )*dΩ ))
end

function solve_crStokes(partition, cells, u_exact, p_exact)
  f(x) = - Δ(u_exact)(x) + ∇(p_exact)(x)
  model = simplexify(CartesianDiscreteModel(partition, cells))

  reffe_u = CrouzeixRaviartRefFE(VectorValue{2,Float64},TRI,1)
  reffe_p = ReferenceFE(lagrangian, Float64, 0; space=:P)
  V = FESpace(model,reffe_u,dirichlet_tags="boundary")
  Q = FESpace(model,reffe_p,conformity=:L2,constraint=:zeromean)
  Y = MultiFieldFESpace([V,Q])

  U = TrialFESpace(V, u_exact)
  P = TrialFESpace(Q)
  X = MultiFieldFESpace([U,P])

  Ω  = Triangulation(model)
  dΩ = Measure(Ω,3)

  a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ
  l((v,q)) = ∫( v⋅f )dΩ

  op = AffineFEOperator(a,l,X,Y)
  uh, ph = solve(op)
  eu = uh-u_exact
  ep = ph-p_exact

  return sqrt(sum( ∫( eu⋅eu )*dΩ )), sqrt(sum( ∫( ep⋅ep )*dΩ ))
end

function conv_test_Stokes(partition,ns,u,p)
    el2u = Float64[]
    el2p = Float64[]
    hs = Float64[]
    for n in ns
        l2u, l2p = solve_crStokes(partition,(n,n),u,p)
        h = 1.0/n
        push!(el2u,l2u)
        push!(el2p,l2p)
        push!(hs,h)
    end
    #println(el2u)
    #println(el2p)
    el2u, el2p, hs
end

function conv_test_Poisson(partition,ns,u)
  el2 = Float64[]
  hs = Float64[]
  for n in ns
  l2 = solve_crScalarPoisson(partition,(n,n),u)
  #println(l2)
  h = 1.0/n
  push!(el2,l2)
  push!(hs,h)
  end
  #println(el2)
  el2, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end


partition = (0,1,0,1)

# Stokes
u_exact(x) = VectorValue( [sin(pi*x[1])^2*sin(pi*x[2])*cos(pi*x[2]), -sin(pi*x[2])^2*sin(pi*x[1])*cos(pi*x[1])] )
p_exact(x) = sin(2*pi*x[1])*sin(2*pi*x[2])

# Scalar Poisson
u_exact_p(x) = sin(2*π*x[1])*sin(2*π*x[2])

# Vector Poisson
# u_exact(x) = VectorValue( [sin(2*π*x[1])*sin(2*π*x[2]), x[1]*(x[1]-1)*x[2]*(x[2]-1)] )

ns = [4,8,16,32,64]

el, hs = conv_test_Poisson(partition,ns,u_exact_p)
# println("Slope L2-norm u_Poisson: $(slope(hs,el))")

elu, elp, hs = conv_test_Stokes(partition,ns,u_exact,p_exact)
@test slope(hs,elu) >= 1.9
@test slope(hs,elp) >= 0.9
@test slope(hs,el)  >= 1.9

#println("Slope L2-norm u_Stokes: $(slope(hs,elu))")
#println("Slope L2-norm p_Stokes: $(slope(hs,elp))")
#println("Slope L2-norm u_Poisson: $(slope(hs,el))")


# for fun to test _get_dfaces_measure in 4D
SIMPL4 = ExtrusionPolytope(tfill(TET_AXIS,Val(4)))
reffe = ReferenceFE(SIMPL4, crouzeix_raviart, 1)
reffec = CrouzeixRaviartRefFE(Float64, SIMPL4, 1)
@test reffe == reffec
@test Conformity(reffe, :L2) == L2Conformity()
@test_throws ErrorException Conformity(reffe, :H1)

end # module

