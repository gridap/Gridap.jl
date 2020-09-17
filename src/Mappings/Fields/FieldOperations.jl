struct FieldOperation{T} <: Mapping
  op::T
end

evaluate!(cache,k::typeof{+},args::NewField...) = GenericField(Operation(k,args...))
evaluate!(cache,k::typeof{-},f::NewField) = GenericField(Operation(k,f))
evaluate!(cache,k::typeof{-},f::NewField,g::NewField) = GenericField(Operation(k,f,g))
evaluate!(cache,k::typeof{inner},f::NewField,g::NewField) = GenericField(Operation(k,f,g))
# add here other meaningful field operations

function evaluate!(cache,k::FieldOperation,args::NewField...)
  GenericField(Operation(k.op,args...))
end

struct Operation{T,F} <: Mapping
  op::T
  args::F
  function Operation(op::Function,args...)
    new{typeof(op),typeof(args)}(op,args)
  end
end

function return_cache(f::Operation,x)
  cargs = return_caches(f.args,x)
  cop = return_cache(op,evaluate!(cargs,f.args,x))
  (cop, cargs)
end

function evaluate!(cache,f::Operation,x)
  cop, cargs = cache
  argsx = evaluate!(cargs,f.args,x)
  evaluate!(cop,f.arg,x)
  # broadcast(f.op,argsx...) # TODO use cache
end

# Implement product rules e.g.
function gradient(a::Operation{typeof(+)})
  Operation(+,map(gradient,a.args))
end

# Operations

# @santiagobadia : My main concern is here, there is no way to decide whether
# the operation is a Field or FieldArray, it depends on the operator, so
# I would be more specific ...
# TODO operations between Fields and arrays of fields
function evaluate!(cache,k::FieldOperation,args::FieldOrFieldArray...)
  GenericFieldArray(Operation(k.op,args...))
end

# TODO this is the difficult part
Base.size(a::Operation) = @notimplemented
Base.getindex(a::Operation,i::Integer...) = @notimplemented
Base.ndims(a::Operation) = @notimplemented
Base.eltype(a::Operation) = @notimplemented
Base.IndexStyle(::Type{<:Operation}) = @notimplemented

# for op in (:+,:-,:tr, :transpose, :adjoint, :symmetric_part)
#   @eval begin

#     function ($op)(f::Field)
#       operate_fields($op,f)
#     end

#     #function ($op)(f::AbstractArray{<:Field})
#     #  operate_arrays_of_fields($op,f)
#     #end

#   end
# end

# for op in (:+,:-,:*,:inner,:outer,:dot)
#   @eval begin

#     function ($op)(f::Field,g::Field)
#       operate_fields($op,f,g)
#     end

#     #function ($op)(f::AbstractArray{<:Field},g::AbstractArray{<:Field})
#     #  operate_arrays_of_fields($op,f,g)
#     #end

#   end
# end

# function operate_fields(op::Function,args...)
#   k = FieldOpKernel(op)
#   apply_kernel_to_field(k,args...)
# end
