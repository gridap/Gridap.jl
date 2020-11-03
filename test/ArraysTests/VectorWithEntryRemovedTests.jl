module VectorWithEntryRemovedTests

using Test
using Gridap.Arrays

a = collect(1:10)
for i=1:length(a)
   b    = VectorWithEntryRemoved(a,i)
   test_array(vcat(a[1:i-1],a[i+1:end]),b)
end

for i=1:length(a)
  b    = VectorWithEntryRemoved(a,i)
  b[rand(1:length(a)-1)]=rand(Int)
  test_array(vcat(a[1:i-1],a[i+1:end]),b)
end

end # module
