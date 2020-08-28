struct FunctionMapping{F} <: NewMapping
  f::F
end

evaluate!(cache,k::FunctionMapping,x...) = k.f(x...)

function return_type(k::FunctionMapping,x...)
  Ts = map(typeof,x)
  return_type(k.f,Ts...)
end

# @santiagobadia : Why not just this?
#function return_type(k::FunctionMapping,x...)
#  typeof(evaluate(k,x...))
#end
