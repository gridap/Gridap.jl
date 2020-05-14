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

add_entries!(a,b)
@test all( a .== ( c .+ b) )

a = rand(6)
c = copy(a)
b = rand(6)

add_entries!(a,b,-)
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

end # module
