
"""
    set_execution_mode(new_mode::String)

Sets the execution mode to either "debug" or "performance", which controls 
the behavior of the @check macro within the Gridap package.
"""
function set_execution_mode(new_mode::String)
  if !(new_mode in ("debug", "performance"))
      throw(ArgumentError("Invalid execution mode: \"$(new_mode)\""))
  end

  # Set it in our runtime values, as well as saving it to disk
  @set_preferences!("execution_mode" => new_mode)
  @info("New execution mode set; restart your Julia session for this change to take effect!")
end

@static if VERSION >= v"1.6"
  const execution_mode = @load_preference("execution_mode", "debug")
else
  const execution_mode = "debug"
end
