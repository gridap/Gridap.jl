module HelperFunctionsTests

using Test
using Gridap.Helpers

# tfill

for D = 0:5
  @test tfill(2,Val(D)) == tuple(fill(2,D)...)
end

# get_val_parameter

@test 1 == get_val_parameter(Val{1}())

# findfirstvalue

v = rand(1:20, 10)
for i in 1:20
  equals_i = ==(i)
  i_or_nothing = i ∈ v ? i : nothing
  @test findfirstvalue(equals_i, v) == i_or_nothing
end

# public_names_in_md

module MockModule
  const C1 = nothing
  const C2 = C1 # Alias
  const C3 = nothing

  export C1
  export C2
  #public C3  # would only work if we stop testing on VERSION < 1.11...
end

change_link=Dict(:C2 => "C1")
s = public_names_in_md(MockModule; change_link)
@test s == """
  ### Exported names

   [`C1`](@ref), [`C2`](@ref C1),"""

  # ### Other public names

  #  [`C3`](@ref),
  # """

end # module
