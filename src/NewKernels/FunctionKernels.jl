struct FunctionKernel{F} <: NewKernel
  f::F
end

evaluate!(cache,k::FunctionKernel,x...) = k.f(x...)

function return_type(k::FunctionKernel,x...)
  Ts = map(typeof,x)
  return_type(k.f,Ts...)
end

# @santiagobadia : Why not just this?
#function return_type(k::FunctionKernel,x...)
#  typeof(evaluate(k,x...))
#end
