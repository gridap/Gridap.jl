
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

@generated function viewtosize(a::Array{T,N},s::NTuple{N,Int}) where {T,N}
    @assert N > 0
    str = join([ ", 1:s[$i]" for i in 1:N ])
    Meta.parse("view(a$str)")
end

# Unicode aliases
const Base.:∘(f,g) = compose(f,g)

const ∇(f) = gradient(f)
