# With Robin boundary conditions

# Manufactured Robin function
rfun(x) = 1.0 + ufun(x)

# Construct the FEspace
diritags = [1,2,3,4,5,7]
fespace = H1ConformingFESpace(Float64,model,order,diritags)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Setup integration on Robin boundary
robintags = [6,]
rtrian = BoundaryTriangulation(model,robintags)
rquad = CellQuadrature(rtrian,degree=2)

# Object describing Robin function
rfield = CellField(rtrian,rfun)

# Integrands of the Robin BC
r(v) = inner(v,rfield)
m(v,u) = inner(v,u)

# Define Robin terms
t_ΓR = AffineFETerm(m,r,rtrian,rquad)

# Dummy term that includes jumps.
# Only for testing purposes since the shape
# functions are continuous in this example.
# Thus, including this term should not change the solution
tags = [9,]
strian = SkeletonTriangulation(model,tags)
squad = CellQuadrature(strian,degree=2)
j(v,u) = inner(jump(v),jump(u))
t_ΓS = LinearFETerm(j,strian,squad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define FE problem
op = LinearFEOperator(V,U,assem,t_Ω,t_ΓN,t_ΓR,t_ΓS)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
e = u - uh

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])

@test el2 < 1.e-8
@test eh1 < 1.e-8
