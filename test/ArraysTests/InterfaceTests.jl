module InterfaceTests

using Test
using Gridap.Arrays

a = rand(20,12)

test_array(a,a)
test_array(a,a,≈)

# @test array_caches() == ()
# @test getitems!((),(),1) == ()

t = (1,2.0,1im,[3,4])
@test testvalue(t) == (0,0.0,0im,Int[])

for n in 0:10
  t = ntuple(i->rand(),n)
  @test testvalue(t) == ntuple(i->0.0,n)
end

end # module
