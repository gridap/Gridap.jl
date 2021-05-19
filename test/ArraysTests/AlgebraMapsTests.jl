module AlgebraMapsTests

using Test
using Gridap.Arrays

a = rand(3,4)
b = rand(4,5)

c = a*b
test_map(c,*,a,b)
cache = return_cache(*,a,b)
@test evaluate!(cache,*,a,b) === evaluate!(cache,*,a,b)

a = rand(3,4)
b = rand(4)

c = a*b
test_map(c,*,a,b)
cache = return_cache(*,a,b)
@test evaluate!(cache,*,a,b) === evaluate!(cache,*,a,b)

a = 3
b = 4

c = a*b
test_map(c,*,a,b)

k! = AddEntriesMap(+)

A = zeros(4,5)
k!(A,fill(3.0,2,2),[1,3],[2,1])
k!(A,fill(3.0,2,2),[1,3],[2,1])
@test A[[1,3],[2,1]] == fill(6.0,2,2)

A = zeros(4)
k!(A,fill(3.0,2),[1,3])
k!(A,fill(3.0,2),[1,3])
@test A[[1,3]] == fill(6.0,2)

k! = TouchEntriesMap()
A = zeros(4,5)
k!(A,fill(3.0,2,2),[1,1],[2,1])
@test A == zeros(4,5)

A = zeros(4)
k!(A,fill(3.0,2),[1,1])
@test A == zeros(4)

end # module
