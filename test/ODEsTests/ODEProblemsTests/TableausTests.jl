module TableausTests

using Test

using Gridap
using Gridap.ODEs
matrix = [
  1 2 3
  4 5 6
  7 8 9
]
@test ODEs._tableau_type(matrix) == FullyImplicitTableau

matrix = [
  1 0 3
  4 5 0
  7 8 0
]
@test ODEs._tableau_type(matrix) == FullyImplicitTableau

matrix = [
  1 0 0
  4 5 0
  7 8 0
]
@test ODEs._tableau_type(matrix) == DiagonallyImplicitTableau

matrix = [
  0 0 0
  4 0 0
  7 8 0
]
@test ODEs._tableau_type(matrix) == ExplicitTableau

for tableauname in available_tableaus
  ButcherTableau(tableauname)
  tableau = eval(Meta.parse("ODEs." * string(tableauname) * "()"))
  ButcherTableau(tableau, Float64)
end

end # module TableausTests
