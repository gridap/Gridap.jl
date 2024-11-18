module JLD2Tests

using Test
using Gridap.Helpers
import Gridap.Io: to_jld2_file
import Gridap.Io: from_jld2_file

struct Foo <: GridapType
  bar::Int
end

bar = 3
foo = Foo(bar)

d = mktempdir()
f = joinpath(d,"foo.jld2")

to_jld2_file(foo,f)
@test foo == from_jld2_file(typeof(foo),f)

f = joinpath(d,"dict.jld2")
foo = Dict("a"=>Int32(1),2=>Int(3),4.0=>Float32(5),"six"=>Float64(7),:s=>"Symbol")

to_jld2_file(foo,f)
@test foo == from_jld2_file(typeof(foo),f)

end # module
