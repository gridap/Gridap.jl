module KeyToValMapsTests

using Test
using Gridap.Arrays
using Gridap.ReferenceFEs

using Gridap
using FillArrays

keys = [1,2,3,4,5]
key_to_val(i) = i*ones(10)

m = KeyToValMap(key_to_val)
r = lazy_map(m,keys)

r[1]
c = return_cache(r)

_r = key_to_val.(keys)
@test all(_r .== collect(r))


key_to_val(i) = LagrangianRefFE(Float64,HEX,i)

m = KeyToValMap(key_to_val)
keys = [1,2,3,4,5,4,3,2,1]
r = lazy_map(m,keys)
_r = key_to_val.(keys)
@test all(get_order.(_r)-get_order.(r) .== 0)

end #module

# @time key_to_val(1)
# @time for i in r end
# @time for i in keys key_to_val(i) end

# @time key_to_val(5)
# dict = Dict(keys .=> key_to_val.(keys))
# # fetching from dict is fast
# @time dict[5]
# @time r[5]
# a = getkey(dict,5,1)
# a
# return_cache(m,keys)

# haskey(dict,5)
# get(dict,5,nothing)
# get!(dict,5,nothing)
