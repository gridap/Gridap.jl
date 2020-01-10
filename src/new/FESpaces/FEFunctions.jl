
"""
"""
abstract type FEFunction <: GridapType end

"""
"""
function get_free_values(f::FEFunction)
  @abstractmethod
end

"""
"""
function get_fe_space(f::FEFunction)
  @abstractmethod
end

"""
"""
function test_fe_function(f::FEFunction)
  free_values = get_free_values(f)
  fe_space = get_fe_space(f)
  @test length(free_values) == num_free_dofs(fe_space)
end

