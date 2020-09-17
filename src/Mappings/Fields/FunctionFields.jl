function_field(f::Function) =  GenericField(f)

const FunctionField{F} = GenericField{F<:Function}

function return_cache(f::FunctionField,x::AbstractArray{<:Point})
  nx = length(x)
  Te = eltype(x)
  c = zeros(return_type(f.f,Te),nx)
  CachedArray(c)
end

function evaluate!(c,f::FunctionField,x::AbstractArray{<:Point})
  nx = length(x)
  setsize!(c,(nx,))
  for i in eachindex(x)
    c[i] = f.f(x[i])
  end
  c
end

function gradient(f::FunctionField)
  FunctionField(gradient(f.f))
end
