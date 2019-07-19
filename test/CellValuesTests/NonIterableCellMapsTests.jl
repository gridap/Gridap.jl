module NonIterableCellMapsTests

using Test
using Gridap

struct Foo <: NonIterableCellMap{Point{2},1,Int,1} end

foo = Foo()

@test isa(foo,CellMap)

try
  iterate(foo)
catch LoadError
end

try
  iterate(foo,2)
catch LoadError
end

end # module
