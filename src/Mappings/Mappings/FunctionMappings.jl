struct FunctionMapping{F} <: Mapping
  f::F
end

evaluate!(cache,k::FunctionMapping,x...) = k.f(x...)

function return_type(k::FunctionMapping,x...)
  Ts = map(typeof,x)
  return_type(k.f,Ts...)
end

evaluate(f::Function,x...) = evaluate(FunctionMapping(f),x...)
evaluate(T,f::Function,x...) = evaluate(T,FunctionMapping(f),x...)

# @santiagobadia : Why not just this?
#function return_type(k::FunctionMapping,x...)
#  typeof(evaluate(k,x...))
#end
