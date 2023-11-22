module TableausTests

using Test

using Gridap
using Gridap.ODEs

matrix = [
  1 2 3
  4 5 6
  7 8 9
]
@test ODEs._butcher_tableau_type(matrix) == ODEs.FullyImplicitTableau

matrix = [
  1 0 0
  4 5 0
  7 8 0
]
@test ODEs._butcher_tableau_type(matrix) == ODEs.DiagonallyImplicitTableau

matrix = [
  0 0 0
  4 0 0
  7 8 0
]
@test ODEs._butcher_tableau_type(matrix) == ODEs.ExplicitTableau

for tableauname in available_tableaus
  ButcherTableau(tableauname)
  ButcherTableau(eval(Meta.parse("ODEs." * string(tableauname) * "()")), Float64)
end

end # module TableausTests
