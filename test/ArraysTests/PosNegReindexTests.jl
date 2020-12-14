module PosNegReindexTests

using Test
using Gridap.Arrays
using Gridap.TensorValues

for indices in ([1,3,-1,2,-2], PosNegPartition([1,4,2],5))

  values_pos = Float64[40,30,10]
  values_neg = -Float64[40,30]
  c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
  r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices ]
  test_array(c,r)

  subindices = [1,4,2]
  d = lazy_map(Reindex(c),subindices)
  @test d === values_pos

  subindices = [3,5]
  d = lazy_map(Reindex(c),subindices)
  @test d === values_neg

  subindices = [4,2]
  d = lazy_map(Reindex(c),subindices)
  r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
  test_array(d,r)

  subindices = [5,]
  d = lazy_map(Reindex(c),subindices)
  r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
  test_array(d,r)

  subindices = [1,5,3,4]
  d = lazy_map(Reindex(c),subindices)
  r = [  i>0 ? values_pos[i] : values_neg[-i] for i in indices[subindices] ]
  test_array(d,r)
  #print_op_tree(d)

end

# Combination with Broadcasting

p = 1
data = [2,3,-1,3,4,4,3,2,5,4,-3,4]
ptrs = [1,4,4,7,13]
gid_to_val_pos = [p, 2*p, 3*p, -p, p]
gid_to_val_neg = [3*p, -p, p]
lid_to_gid = Table(data,ptrs)
ca = lazy_map(Broadcasting(PosNegReindex(gid_to_val_pos,gid_to_val_neg)),lid_to_gid)
a = Vector{Int}[[2, 3, 3], [], [3, -1, -1], [3, 2, 1, -1, 1, -1]]
test_array( ca, a )

i_to_v=[[1, 2], [5, 6], [1, 5], [2, 6], [2, 3], [6, 7], [3, 7], [3, 4], [7, 8], [4, 8], [9, 10], [5, 9], [6, 10], [10, 11], [7, 11], [11, 12], [8, 12], [13, 14], [9, 13], [10, 14], [14, 15], [11, 15], [15, 16], [12, 16]]
j_to_i=Int[]
lid_to_gid=lazy_map(Reindex(i_to_v),j_to_i)
gid_to_val=VectorValue{2,Float64}[(0.0, 0.25), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (0.0, 0.5), (0.25, 0.5), (0.5, 0.5), (0.75, 0.5), (0.0, 0.75), (0.25, 0.75), (0.5, 0.75), (0.75, 0.75), (0.0, 1.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0)]
f2=lazy_map(Broadcasting(PosNegReindex(gid_to_val,gid_to_val)),lid_to_gid)

# 0-length cases

test_array(f2,Vector{VectorValue{2,Float64}}[])

indices = Int[]
values_pos = Float64[]
values_neg = Float64[]
c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
test_array(c,Float64[])

# If we work with PosNegReindex return also PosNegReindex

for indices in ([1,3,-2,2,-1], PosNegPartition([1,4,2],5), lazy_map(PosNegReindex([1,2,3],[-1,-2]),[1,3,-2,2,-1]))

  a_pos = Float64[40,30,10]
  a_neg = -Float64[40,30]
  a = lazy_map(PosNegReindex(a_pos,a_neg),indices)

  b_pos = Float64[43,50,60]
  b_neg = -Float64[41,30]
  b = lazy_map(PosNegReindex(b_pos,b_neg),indices)

  f_pos = fill(+,3)
  f_neg = fill(-,2)
  f = lazy_map(PosNegReindex(f_pos,f_neg),Function,indices)
  c = lazy_map(evaluate,f,a,b)
  #print_op_tree(c)
  test_array(c,collect(c))
  c_pos = c.maps.value.values_pos
  c_neg = c.maps.value.values_neg
  test_array(c_pos,a_pos+b_pos)
  test_array(c_neg,a_neg-b_neg)

  a_pos = Float64[40,30,10]
  a_neg = -Int[40,30]
  a = lazy_map(PosNegReindex(a_pos,a_neg),Number,indices)
  fun(i) = Float64(i) + 3
  c = lazy_map(fun,a)
  #print_op_tree(c)
  test_array(c,collect(c))
  c_pos = c.maps.value.values_pos
  c_neg = c.maps.value.values_neg
  test_array(c_pos,fun.(a_pos))
  test_array(c_neg,fun.(a_neg))

end

a_i = PosNegPartition([1,],4)
a_pos = Float64[40]
a_neg = -Float64[40,30,20]
a = lazy_map(PosNegReindex(a_pos,a_neg),a_i)

b_i = PosNegPartition([3,2,4],4)
b_pos = Float64[50,60,30]
b_neg = -Float64[41]
b = lazy_map(PosNegReindex(b_pos,b_neg),b_i)

c = lazy_map(+,a,b)
test_array(c,map(+,a,b))



## Testing some cases where PosNegReindex can be type-instable
#
#indices = [1,3,-2,2,-1]
#values_pos = Int[40,30,10]
#values_neg = -Float64[40,30]
#c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
#d = lazy_map(i->Float64(i),c)
#display(d)
#
#d = lazy_map(*,c,ones(Float64,size(c)))
#display(d)
#
#indices = [1,1,-1,1,-1]
#values_pos = [+]
#values_neg = [-]
#c = lazy_map(PosNegReindex(values_pos,values_neg),indices)
#d = lazy_map(evaluate, c, rand(size(c)), rand(size(c)))
#display(d)

end # module
