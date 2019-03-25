
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

flatten(a::Array) = reshape(a,(length(a),))

# Unicode aliases
const Base.:∘(f,g) = compose(f,g)

const ∇(f) = gradient(f)
