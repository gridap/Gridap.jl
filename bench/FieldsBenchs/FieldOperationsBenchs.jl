module FieldOperationsBenchs

using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: MockField, MockBasis
using Gridap.TensorValues
using FillArrays

@inline function loop(a,cache)
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
  end
end

function bench1(n)
  p1 = Point(2,2)
  p2 = Point(4,2)
  p3 = Point(1,3)
  p4 = Point(5,2)
  x = [p1,p2,p3,p4]
  np = length(x)
  va = 3.0
  vb = 3.3
  d = 2
  fa = MockField{d}(va)
  fb = MockField{d}(vb)
  wa = 4.5
  wb = 4.5
  ndofa = 8
  ndofb = 5
  ba = MockBasis{d}(wa,ndofa)
  bb = MockBasis{d}(wb,ndofb)
  l = n
  ax = fill(x,l)
  afa = Fill(fa,l)
  afb = Fill(fb,l)
  aba = Fill(ba,l)
  abb = Fill(bb,l)

  ag_ff = operate_arrays_of_fields(inner,afa,afb)
  agx_ff = evaluate(ag_ff,ax)
  cagx_ff = array_cache(agx_ff)
  @time loop(agx_ff,cagx_ff)

  ag_bf = operate_arrays_of_fields(inner,aba,afb)
  agx_bf = evaluate(ag_bf,ax)
  cagx_bf = array_cache(agx_bf)
  @time loop(agx_bf,cagx_bf)

  #ag_bb = operate_arrays_of_fields(inner,aba,abb)
  #agx_bb = evaluate(ag_bb,ax)
  #cagx_bb = array_cache(agx_bb)
  #@time loop(agx_bb,cagx_bb)

end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ running suite for n = $($n) +++")
    bench1($n)
  end
end


end # module
