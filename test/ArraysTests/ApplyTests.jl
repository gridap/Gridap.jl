module ApplyTests

using Test
using Gridap.Arrays
using FillArrays

a = rand(3,2,4)
test_array(a,a)

a = rand(3,2,4)
a = CartesianIndices(a)
test_array(a,a)

a = rand(3,2)
a = CartesianIndices(a)
c = apply(-,a)
test_array(c,-a)

a = rand(12)
c = apply(-,a)
test_array(c,-a)

a = rand(12)
b = rand(12)
c = apply(-,a,b)
test_array(c,a.-b)

c = apply(Float64,-,a,b)
test_array(c,a.-b)

a = rand(0)
b = rand(0)
c = apply(-,a,b)
test_array(c,a.-b)

a = fill(rand(2,3),12)
b = rand(12)
c = apply(bcast(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),0)
b = rand(0)
c = apply(bcast(-),a,b)
test_array(c,[ai.-bi for (ai,bi) in zip(a,b)])

a = fill(rand(2,3),12)
b = rand(12)
c = apply(bcast(-),a,b)
d = apply(bcast(+),a,c)
e = apply(bcast(*),d,c)
test_array(e,[((ai.-bi).+ai).*(ai.-bi) for (ai,bi) in zip(a,b)])

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
c = array_caches(a,b)
i = 1
v = getitems!(c,(a,b),i)
@test c == (nothing,nothing)
@test v == (a[i],b[i])

a = fill(rand(Int,2,3),12)
b = fill(rand(Int,1,3),12)
ai = testitem(a)
@test ai == a[1]
ai, bi = testitems(a,b)
@test ai == a[1]
@test bi == b[1]

a = fill(rand(Int,2,3),0)
b = fill(1,0)
ai = testitem(a)
@test ai == Array{Int,2}(undef,0,0)
ai, bi = testitems(a,b)
@test ai == Array{Int,2}(undef,0,0)
@test bi == zero(Int)

a = fill(+,10)
x = rand(10)
y = rand(10)
v = apply(a,x,y)
r = [(xi+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)
v = apply(Float64,a,x,y)
test_array(v,r)

a = Fill(bcast(+),10)
x = [rand(2,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

a = Fill(bcast(+),10)
x = [rand(mod(i-1,3)+1,3) for i in 1:10]
y = [rand(1,3) for i in 1:10]
v = apply(a,x,y)
r = [(xi.+yi) for (xi,yi) in zip(x,y)]
test_array(v,r)

# Test the intermediate results caching mechanism

struct ArrayWithCounter{T,N,A,C} <: AbstractArray{T,N}
  array::A
  counter::C
  function ArrayWithCounter(a::AbstractArray{T,N}) where {T,N}
    c = zeros(Int,size(a))
    c[:] .= 0
    new{T,N,typeof(a),typeof(c)}(a,c)
  end
end

Base.size(a::ArrayWithCounter) = size(a.array)

function Base.getindex(a::ArrayWithCounter,i::Integer...)
  a.counter[i...] += 1
  a.array[i...]
end

Base.IndexStyle(::Type{<:ArrayWithCounter{T,N,A}}) where {T,A,N} = IndexStyle(A)

function reset_counter!(a::ArrayWithCounter)
  a.counter[:] .= 0
end

a = ArrayWithCounter(fill(rand(2,3),12))
b = ArrayWithCounter(rand(12))
c = apply(bcast(-),a,b)
d = apply(bcast(+),a,c)
e = apply(bcast(*),d,c)
r = [ (ai.-bi).*(ai.+(ai.-bi)) for (ai,bi) in zip(a,b)]
cache = array_cache(e)
reset_counter!(a)
reset_counter!(b)
for i in 1:length(e)
  ei = getindex!(cache,e,i)
  ei = getindex!(cache,e,i)
  ei = getindex!(cache,e,i)
end

@test all(a.counter .== 2) 
@test all(b.counter .== 1)

l = 10
ai = 1.0
bi = 2.0
a = Fill(ai,l)
b = Fill(bi,l)
c = apply(+,a,b)
r = map(+,a,b)
test_array(c,r)
@test isa(c,Fill)

#f = rand(10)
#a = rand(10)
#@test apply(f,a) === f
#
#f = fill(rand(4),10)
#a = rand(10)
#@test apply(f,a) === f
#
#f = fill(rand(4),10)
#g = fill(rand(4),10)
#a = rand(10)
#h = apply_all((f,g),a)
#@test h[1] === f
#@test h[2] === g

#l = 10
#ai = 1.0
#bi = 2.0
#a = fill(ai,l)
#b = fill(bi,l)
#c = apply(+,a,b)
#d = apply(c,a,b)
#@test c == d
#@test typeof(d) == typeof(c)

end # module
