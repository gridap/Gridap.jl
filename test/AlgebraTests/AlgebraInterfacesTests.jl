module AlgebraInterfacesTests

using Gridap.Algebra
using Test

a = allocate_vector(Vector{Int},10)
@test isa(a,Vector{Int})
@test length(a) == 10

a = allocate_vector(Vector{Int},1:10)
@test isa(a,Vector{Int})
@test length(a) == 10

a = zeros(4,6)

b = allocate_in_range(Vector{Int},a)
@test length(b) == 4

b = allocate_in_domain(Vector{Int},a)
@test length(b) == 6

fill_entries!(b,4)
@test all( b .== 4 )

a = rand(6)
b = rand(6)

copy_entries!(a,b)
@test all( a .== b )

a = rand(6)
c = copy(a)
b = rand(6)

a .+= b
@test all( a .== ( c .+ b) )

a = rand(6)
c = copy(a)
b = rand(6)

a .-= b
@test all( a .== ( c .- b) )

a = rand(6)
c = copy(a)
scale_entries!(a,10)
@test all( a .== 10*c)

a = rand(4,6)
c = rand(4)
d = copy(c)
b = rand(6)

muladd!(c,a,b)

@test all( c .≈ (d .+ a*b ) )

A = Vector{Float64}
n = 10
rows = Base.OneTo(n)
a = nz_counter(A,(rows,))
@test LoopStyle(a) == DoNotLoop()
add_entry!(a,1.0,1)
add_entries!(a,[1.0,-1.0],[1,1])
add_entries!(a,nothing,[1,1])
b = nz_allocation(a)
@test LoopStyle(b) == DoNotLoop()
@test isa(b,Vector{Float64})
@test length(b) == n
add_entry!(b,1.0,1)
add_entry!(b,1.0,1)
add_entry!(b,1.0,4)
add_entries!(b,[1.0,-1.0],[1,1])
add_entries!(b,nothing,[1,1])
r = zeros(n)
r[1] = 2
r[4] = 1
@test b == r
c = create_from_nz(b)
@test c === b
add_entries!(c,[1.0,-1.0],[1,1])
add_entries!(c,nothing,[1,1])

using SparseArrays
using SparseMatricesCSR
for A in (
  SparseMatrixCSC{Float64,Int},
  SparseMatrixCSR{1,Float64,Int},
  SparseMatrixCSR{0,Float64,Int})

  m = 6
  n = 9
  rows = Base.OneTo(m)
  cols = Base.OneTo(n)
  a = nz_counter(A,(rows,cols))
  @test LoopStyle(a) == Loop()
  add_entry!(a,1.0,1,1)
  add_entry!(a,nothing,1,1)
  add_entry!(a,nothing,3,1)
  add_entry!(a,3.0,4,9)
  @test a.nnz == 4
  b = nz_allocation(a)
  @test LoopStyle(b) == Loop()
  @test length(b.I) == a.nnz
  @test length(b.J) == a.nnz
  @test length(b.V) == a.nnz
  add_entry!(b,1.0,1,1)
  add_entry!(b,nothing,1,1)
  add_entry!(b,nothing,3,1)
  add_entry!(b,3.0,4,9)
  c = create_from_nz(b)
  @test LoopStyle(c) == DoNotLoop()
  @test isa(c,A)
  I,J,V = findnz(c)
  @test I == [1,3,4]
  @test J == [1,1,9]
  @test V == Float64[1,0,3]
  add_entry!(c,1.0,1,1)
  add_entry!(c,nothing,1,1)
  add_entry!(c,nothing,3,1)
  add_entry!(c,3.0,4,9)
  
  a = nz_counter(A,(rows,cols))
  add_entries!(a,[1.0 -1.0; -1.0 1.0],[1,-1],[-1,1])
  add_entries!(a,nothing,[1,1],[1,-1])
  @test a.nnz == 3
  b = nz_allocation(a)
  add_entries!(b,[1.0 -1.0; -1.0 1.0],[1,-1],[-1,1])
  add_entries!(b,nothing,[1,1],[1,-1])
  c = create_from_nz(b)
  add_entries!(c,[1.0 -1.0; -1.0 1.0],[1,-1],[-1,1])
  add_entries!(c,nothing,[1,1],[1,-1])
end

for A in (
  SymSparseMatrixCSR{1,Float64,Int},
  SymSparseMatrixCSR{0,Float64,Int})

  m = 6
  n = 6
  rows = Base.OneTo(m)
  cols = Base.OneTo(n)
  a = nz_counter(A,(rows,cols))
  @test LoopStyle(a) == Loop()
  add_entry!(a,1.0,1,1)
  add_entry!(a,nothing,1,1)
  add_entry!(a,nothing,3,1)
  add_entry!(a,nothing,1,3)
  add_entry!(a,3.0,2,2)
  add_entry!(a,3.0,2,6)
  add_entry!(a,3.0,6,2)
  @test a.nnz == 5
  b = nz_allocation(a)
  add_entry!(b,1.0,1,1)
  add_entry!(b,nothing,1,1)
  add_entry!(b,nothing,3,1)
  add_entry!(b,nothing,1,3)
  add_entry!(b,3.0,2,2)
  add_entry!(b,3.0,2,6)
  add_entry!(b,3.0,6,2)
  @test LoopStyle(b) == Loop()
  @test length(b.I) == a.nnz
  @test length(b.J) == a.nnz
  @test length(b.V) == a.nnz
  c = create_from_nz(b)
  @test LoopStyle(c) == DoNotLoop()
  @test isa(c,A)
  I,J,V = findnz(c)
  @test I == [1,1,2,2,3,4,5,6]
  @test J == [1,3,2,6,3,4,5,6]
  @test V == Float64[1,0,3,3,0,0,0,0]
  add_entry!(c,1.0,1,1)
  add_entry!(c,nothing,1,1)
  add_entry!(c,nothing,3,1)
  add_entry!(c,nothing,1,3)
  add_entry!(c,3.0,2,2)
  add_entry!(c,3.0,2,6)
  add_entry!(c,3.0,6,2)
  
  a = nz_counter(A,(rows,cols))
  add_entries!(a,[1.0 -1.0; -1.0 1.0],[1,-1],[-1,1])
  add_entries!(a,nothing,[1,1],[1,-1])
  @test a.nnz == 3
  b = nz_allocation(a)
  add_entries!(b,[1.0 -1.0; -1.0 1.0],[1,-1],[-1,1])
  add_entries!(b,nothing,[1,1],[1,-1])
  c = create_from_nz(b)
  add_entries!(c,[1.0 -1.0; -1.0 1.0],[1,-1],[-1,1])
  add_entries!(c,nothing,[1,1],[1,-1])
end

end # module
