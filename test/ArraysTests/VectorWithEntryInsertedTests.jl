module VectorWithEntryInsertedTests

using Test
using Gridap.Arrays

a = collect(1:10)
for i=1:length(a)
   b = VectorWithEntryInserted(a,i,i)
   test_array(vcat(a[1:i],[i],a[i+1:end]),b)
end

for i=1:length(a)
  b = VectorWithEntryInserted(a,i,i)
  pos = rand(1:length(b))
  val = rand(Int)
  b[pos]=val
  if (pos == i)
     test_array(vcat(a[1:i-1],[val],a[i:end]),b)
  else
     test_array(vcat(a[1:i-1],[i],a[i:end]),b)
  end
end

end # module
