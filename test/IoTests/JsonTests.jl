module JsonTests

using Test
using Gridap.Helpers
using Gridap.Io
import Gridap.Io: to_dict
import Gridap.Io: from_dict

using JSON

struct Foo <: GridapType
  bar::Int
end

function to_dict(foo::Foo)
  dict = Dict{Symbol,Any}()
  dict[:bar] = foo.bar
  dict
end

function from_dict(::Type{Foo},dict::Dict{Symbol,Any})
  bar = dict[:bar]
  Foo(bar)
end

bar = 3
foo = Foo(bar)

@test foo == from_dict(Foo,to_dict(foo))
@test foo == from_json(Foo,to_json(foo))

d = mktempdir()
f = joinpath(d,"foo.json")
to_json_file(foo,f)
@test foo == from_json_file(Foo,f)
@test foo == from_json(Foo,JSON.json(foo))

end # module




