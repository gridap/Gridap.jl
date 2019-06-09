module Helpers

export @abstractmethod
export @notimplemented
export @notimplementedif
export @unreachable

macro abstractmethod()
  quote
    error("This function belongs to an interface definition and cannot be used.")
  end
end

macro notimplemented()
  quote
    error("This function in not yet implemented")
  end
end

macro notimplementedif(condition)
  quote
    if $(esc(condition))
      @notimplemented
    end
  end
end

macro unreachable()
  quote
    error("This line of code cannot be reached")
  end
end

end # module Helpers
