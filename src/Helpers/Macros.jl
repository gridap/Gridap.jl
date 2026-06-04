
"""
    @abstractmethod

Macro used in generic functions that must be overloaded by derived types.
"""
macro abstractmethod(message="This method belongs to an abstract interface definition and cannot be used.")
  quote
    error($(esc(message)))
  end
end

"""
    @notimplemented
    @notimplemented "Error message"

Macro used to raise an error, when something is not implemented.
"""
macro notimplemented(message="This method is not yet implemented.")
  quote
    error($(esc(message)))
  end
end

"""
    @notimplementedif condition
    @notimplementedif condition "Error message"

Macro used to raise an error if the `condition` is true
"""
macro notimplementedif(condition,message="This method is not yet implemented.")
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

"""
    @check_inferred expression

Macro used unsure that `expression` is inferrable (passes [`Test.@inferred`](@ref)).
The check gets deactivated when running Gridap in performance mode.
"""
macro check_inferred(ex)
  ex.head != :call && error("@check_inferred requires a call expression, got $(string(ex))")
  ex_str = string(ex)
  args_symbs = map(ex.args[2:end]) do arg_ex
    string(arg_ex)
  end
  args_ex = map(ex.args[2:end]) do arg_ex
    esc(arg_ex)
  end

  @static if execution_mode == "debug"
    quote
      try
        $(esc( :(Helpers.Test.@inferred $(ex)) ))
      catch e
        if e isa InferrabilityException
          # This ensures that we track the innermost inferrability failure of the tree
          rethrow(e)

        elseif e isa ErrorException && contains(e.msg , "does not match inferred return type")
          args_symbs = $args_symbs
          args_types = map(string∘typeof, tuple($(args_ex...)))
          args_vals = tuple($(args_ex...))
          ids = 1:length(args_symbs)
          args_str = map(ids, args_symbs, args_types, args_vals) do i, symb, T, val
            "    $symb::$T = $val\n"
          end

          msg = """
          InferrabilityException: The following call is not type-inferrable:
              $($ex_str)

          Call arguments:
          $(join(args_str, "\n"))
          """
          rethrow(InferrabilityException(msg))

        else
          rethrow()
        end
      end # catch
    end # quote
  else
    esc(ex)
  end
end

struct InferrabilityException <: Exception
  msg
end

Base.showerror(io::IO, e::InferrabilityException) = print(io, e.msg)

