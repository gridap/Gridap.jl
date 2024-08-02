
"""
    @abstractmethod

Macro used in generic functions that must be overloaded by derived types.
"""
macro abstractmethod(message="This function belongs to an interface definition and cannot be used.")
  quote
    error($(esc(message)))
  end
end

"""
    @notimplemented
    @notimplemented "Error message"

Macro used to raise an error, when something is not implemented.
"""
macro notimplemented(message="This function is not yet implemented")
  quote
    error($(esc(message)))
  end
end

"""
    @notimplementedif condition
    @notimplementedif condition "Error message"

Macro used to raise an error if the `condition` is true
"""
macro notimplementedif(condition,message="This function is not yet implemented")
  quote
    if $(esc(condition))
      @notimplemented $(esc(message))
    end
  end
end

"""
    @unreachable
    @unreachable "Error message"

Macro used to make sure that a line of code is never reached.
"""
macro unreachable(message="This line of code cannot be reached")
  quote
    error($(esc(message)))
  end
end

"""
  @check condition
  @check condition "Error message"

Macro used to make sure that condition is fulfilled, like `@assert`
but the check gets deactivated when running Gridap in performance mode.
"""
macro check(test,msg="A check failed")
  @static if execution_mode == "debug"
    quote
      @assert $(esc(test)) $(esc(msg))
    end
  else
    quote
      nothing
    end
  end
end
