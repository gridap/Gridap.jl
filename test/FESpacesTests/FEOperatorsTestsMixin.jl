# Construct the FEspace
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,order=order*2)

# Define forms
a(v,u) = inner(∇(v), ∇(u))
b(v) = inner(v,bfun)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(a,b,V,U,assem,trian,quad)

# Define the FESolver
ls = LUSolver()
solver = LinearFESolver(ls)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
u = CellField(trian,ufun)
e = u - uh

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = a(u,u) + l2(u)

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])

# Further tests

@test TrialFESpace(op) === U
@test TestFESpace(op) === V

zh = zero(U)
res = apply(op,zh)
r2 = similar(res)
apply!(r2,op,zh)
@test res ≈ r2

cache = solve!(zh,solver,op)
@test free_dofs(zh) ≈ free_dofs(uh)

zh = zero(U)
solve!(zh,solver,op,cache)
@test free_dofs(zh) ≈ free_dofs(uh)

# With Neumann BCs

# Manufactured Neumann function
gfun(x) = 1.0

# Construct the FEspace
order = 1
diritags = [1,2,3,4,5,6,7]
fespace = ConformingFESpace(Float64,model,order,diritags)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Setup integration on Neumann boundary
neumanntags = [8,]
btrian = BoundaryTriangulation(model,neumanntags)
bquad = CellQuadrature(btrian,order=order*2)

# Integrand of the Neumann BC
g(v) = inner(v,gfun)

# Define weak form terms
t_Ω = AffineFETerm(a,b,trian,quad)
t_ΓN = FESource(g,btrian,bquad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define FE problem
op = LinearFEOperator(V,U,assem,t_Ω,t_ΓN)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
e = u - uh

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])
