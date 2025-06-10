module HelperFunctionsTests

using Test
using Gridap.Helpers

for D = 0:5
  @test tfill(2,Val(D)) == tuple(fill(2,D)...)
end

@test 1 == get_val_parameter(Val{1}())

module MockModule
  const C1 = nothing
  const C2 = C1 # Alias
  const C3 = nothing

  export C1
  export C2
  public C3
end

change_link=Dict(:C2 => "C1")
s = public_names_in_md(MockModule; change_link)
@test s == """
  ### Exported names

   [`C1`](@ref), [`C2`](@ref C1),

  ### Other public names

   [`C3`](@ref),
  """

end # module
