
"""
    set_execution_mode(new_mode::String)

Sets the execution mode to either "debug" or "performance", which controls 
the behavior of the @check macro within the Gridap package.

- Debug mode (default): The @check macro will be active, which activates consistency checks 
  within the library. This mode is recommended for development and debugging purposes.

- Performance mode: The @check macro will be deactivated. This mode is recommended for 
  production runs, where no errors are expected.

Pre-defined functions [`set_debug_mode`](@ref) and [`set_performance_mode`](@ref) are also available.
Feature only available in Julia 1.6 and later due to restrictions from `Preferences.jl`.
"""
function set_execution_mode(new_mode::String)
  if !(new_mode in ("debug", "performance"))
    throw(ArgumentError("Invalid execution mode: \"$(new_mode)\""))
  end

  # Set it in our runtime values, as well as saving it to disk
  @set_preferences!("execution_mode" => new_mode)
  @info("New execution mode set; restart your Julia session for this change to take effect!")
end

"""
    set_debug_mode()

  Equivalent to `set_execution_mode("debug")`.
"""
set_debug_mode() = set_execution_mode("debug")

"""
    set_performance_mode()

  Equivalent to `set_execution_mode("performance")`.
"""
set_performance_mode() = set_execution_mode("performance")

@static if VERSION >= v"1.6"
  const execution_mode = @load_preference("execution_mode", "debug")
else
  const execution_mode = "debug"
end
