using Gridap
import Gridap: ∇

# Define manufactured functions
ufun(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0,0.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = 0.0

# Define forms of the problem
a(v,u) = inner(∇(v), ∇(u))

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = a(u,u) + l2(u)

function main(n)

  # Construct the discrete model
  println("model")
  @time model = CartesianDiscreteModel(partition=(n,n,n))
  
  # Construct the FEspace
  order = 1
  diritag = "boundary"
  println("fespace")
  @time fespace = ConformingFESpace(Float64,model,order,diritag)
  
  # Define test and trial spaces
  println("V")
  @time V = TestFESpace(fespace)

  println("U")
  @time U = TrialFESpace(fespace,ufun)
  
  # Define integration mesh and quadrature
  println("trian")
  @time trian = Triangulation(model)

  println("quad")
  @time quad = CellQuadrature(trian,order=2)
  
  # Define the source term
  println("bfield")
  @time bfield = CellField(trian,bfun)
  
  # Define Assembler
  println("assem")
  @time assem = SparseMatrixAssembler(V,U)
  
  # Define the FEOperator
  b(v) = inner(v,bfield)
  println("op")
  @time op = LinearFEOperator(a,b,V,U,assem,trian,quad)
  
  # Define the FESolver
  ls = LUSolver()
  solver = LinearFESolver(ls)
  
  # Solve!
  println("uh")
  @time uh = solve(solver,op)
  
  # Define exact solution and error
  println("u")
  @time u = CellField(trian,ufun)

  println("e")
  @time e = u - uh
  
  # Compute errors
  println("el2")
  @time el2 = sqrt(sum( integrate(l2(e),trian,quad) ))

  println("eh1")
  @time eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
  
  @assert el2 < 1.e-8
  @assert eh1 < 1.e-8

end

