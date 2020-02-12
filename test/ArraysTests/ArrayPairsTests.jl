module ArrayPairsTests

using Gridap.Arrays

a = collect(1:10)
b = collect(11:20)
r = [(i,i+10) for i in 1:10]

c = pair_arrays(a,b)
test_array(c,r)

end # module
