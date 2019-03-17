
macro abstractmethod()
  quote
    error("This function belongs to an interface definition and cannot be used.")
  end
end

