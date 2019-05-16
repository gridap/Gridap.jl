module Helpers

export @abstractmethod
export @notimplemented
export @notimplementedif
export @unreachable
export viewtosize
export rewind_ptrs!
export length_to_ptrs!

import Gridap: flatten

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

flatten(a::Array) = reshape(a,(length(a),))

@generated function viewtosize(a::Array{T,N},s::NTuple{N,Int}) where {T,N}
    @assert N > 0
    str = join([ ", 1:s[$i]" for i in 1:N ])
    Meta.parse("view(a$str)")
end

function rewind_ptrs!(ptrs::AbstractArray{T,1}) where T
  @inbounds for i in (length(ptrs)-1):-1:1
    ptrs[i+1] = ptrs[i]
  end
  ptrs[1] = 1
end

function length_to_ptrs!(ptrs::AbstractArray{T,1}) where T
  ptrs[1] = 1
  @inbounds for i in 1:(length(ptrs)-1)
    ptrs[i+1] += ptrs[i]
  end
end

end # module Helpers
