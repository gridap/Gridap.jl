struct FunctionField{F} <: Field
  f::F
end

function_field(f::Function) =  FunctionField(f)

function field_cache(f::FunctionField,x)
  nx = length(x)
  Te = eltype(x)
  c = zeros(return_type(f.f,Te),nx)
  CachedArray(c)
end

function evaluate_field!(c,f::FunctionField,x)
  nx = length(x)
  setsize!(c,(nx,))
  for i in eachindex(x)
    c[i] = f.f(x[i])
  end
  c
end

function field_gradient(f::FunctionField)
  FunctionField(gradient(f.f))
end
