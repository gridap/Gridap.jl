module MacrosTests

using Test
using Gridap.Helpers

y = 0
f(x) = y+x
@test_throws Helpers.InferrabilityException @check_inferred f(0)

let y=y
g(x) = y+x
@test_nowarn @check_inferred g(0)
end

const z = 0
h(x) = z+x
@test_nowarn @check_inferred h(0)

end # module
