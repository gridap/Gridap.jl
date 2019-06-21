module CellValuesOperationsTests

using Test
using TensorValues
using Gridap

using Gridap.CachedArrays

using ..CellValuesMocks
using ..MapsMocks

l = 10

# CellNumber Operations

ax = [1.2, VectorValue(1.3,2.0), VectorValue(1.3,2.0)]
bx = [1  , VectorValue(1.1,2.2), 1                   ]

for op in (:+,:-)
  @eval begin

    for a in ax
      v = TestIterCellValue(a,l)
      w = $op(v)
      o = [ $op(vi) for vi in v ] 
      test_iter_cell_value(w,o)
    end

    for a in ax
      v = TestIndexCellValue(a,l)
      w = $op(v)
      o = [ $op(vi) for vi in v ] 
      test_index_cell_value(w,o)
    end

    for (a,b) in zip(ax,bx)
      u = TestIterCellValue(a,l)
      v = TestIterCellValue(b,l)
      w = $op(u,v)
      o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
      test_iter_cell_value(w,o)
    end

    for (a,b) in zip(ax,bx)
      u = TestIndexCellValue(a,l)
      v = TestIndexCellValue(b,l)
      w = $op(u,v)
      o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
      test_index_cell_value(w,o)
    end

  end
end

a = TensorValue(1.0,0.0,1.0,2.0)

for op in (:inv,:det,:meas)
  @eval begin

    v = TestIterCellValue(a,l)
    w = $op(v)
    o = [ $op(vi) for vi in v ] 
    test_iter_cell_value(w,o)

    v = TestIndexCellValue(a,l)
    w = $op(v)
    o = [ $op(vi) for vi in v ] 
    test_index_cell_value(w,o)

  end
end

a = VectorValue(1.0,0.0,1.0)
b = VectorValue(1.0,4.0,-1.0)

for op in (:inner,:outer)
  @eval begin

    u = TestIterCellValue(a,l)
    v = TestIterCellValue(b,l)
    w = $op(u,v)
    o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
    test_iter_cell_value(w,o)

    u = TestIndexCellValue(a,l)
    v = TestIndexCellValue(b,l)
    w = $op(u,v)
    o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
    test_index_cell_value(w,o)

  end
end

a = TensorValue(1.0,0.0,1.0,2.0)
b = VectorValue(1.0,4.0)

for op in (:*,:\)
  @eval begin

    u = TestIterCellValue(a,l)
    v = TestIterCellValue(b,l)
    w = $op(u,v)
    o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
    test_iter_cell_value(w,o)

    u = TestIndexCellValue(a,l)
    v = TestIndexCellValue(b,l)
    w = $op(u,v)
    o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
    test_index_cell_value(w,o)

  end
end

a = TensorValue(1.0,0.0,1.0,2.0)
b = 2

for op in (:/,)
  @eval begin

    u = TestIterCellValue(a,l)
    v = TestIterCellValue(b,l)
    w = $op(u,v)
    o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
    test_iter_cell_value(w,o)

    u = TestIndexCellValue(a,l)
    v = TestIndexCellValue(b,l)
    w = $op(u,v)
    o = [ $op(ui,vi) for (ui,vi) in zip(u,v) ] 
    test_index_cell_value(w,o)

  end
end

u = TestIterCellValue(a,l)
v = TestIterCellValue(a,l)
@test u == v
@test u ≈ v

# CellArray Operations

v1 = VectorValue(2,3)
v2 = VectorValue(3,2)
v3 = VectorValue(1,2)

ax = [rand(1,3,4), 1.0      , [v1,v2,v3], [v1,v2,v3], 1         ]
bx = [rand(2,3,1), rand(2,3), [v2,v3,v1], v1        , [v2,v3,v1]]
cx = [rand(1,3,4), [v1,v2,v3]]

for op in (:+,:-)
  @eval begin

    for a in cx

      v = TestIterCellValue(a,l)
      w = $op(v)
      o = [ CachedArray(broadcast($op,vi)) for vi in v ] 
      test_iter_cell_array(w,o)

      v = TestIndexCellValue(a,l)
      w = $op(v)
      o = [ CachedArray(broadcast($op,vi)) for vi in v ] 
      test_index_cell_array(w,o)

    end

    for (a,b) in zip(ax,bx)

      u = TestIterCellValue(a,l)
      v = TestIterCellValue(b,l)
      o = [ CachedArray(broadcast($op,ui,vi)) for (ui,vi) in zip(u,v) ] 
      w = $op(u,v)
      test_iter_cell_array(w,o)

      u = TestIndexCellValue(a,l)
      v = TestIndexCellValue(b,l)
      o = [ CachedArray(broadcast($op,ui,vi)) for (ui,vi) in zip(u,v) ] 
      w = $op(u,v)
      test_index_cell_array(w,o)

    end

  end
end

ai = TensorValue(1.0,0.0,1.0,2.0)
a = fill(ai,2,3)

for op in (:inv,:det,:meas)
  @eval begin

  v = TestIterCellValue(a,l)
  w = $op(v)
  o = [ CachedArray(broadcast($op,vi)) for vi in v ] 
  test_iter_cell_array(w,o)

  v = TestIndexCellValue(a,l)
  w = $op(v)
  o = [ CachedArray(broadcast($op,vi)) for vi in v ] 
  test_index_cell_array(w,o)

  end
end

ai = VectorValue(1.0,0.0,1.0)
bi = VectorValue(1.0,4.0,-1.0)
a = fill(ai,2,3)
b = fill(bi,1,3)

for op in (:inner,:outer)
  @eval begin

    u = TestIterCellValue(a,l)
    v = TestIterCellValue(b,l)
    o = [ CachedArray(broadcast($op,ui,vi)) for (ui,vi) in zip(u,v) ] 
    w = $op(u,v)
    test_iter_cell_array(w,o)

    u = TestIndexCellValue(a,l)
    v = TestIndexCellValue(b,l)
    o = [ CachedArray(broadcast($op,ui,vi)) for (ui,vi) in zip(u,v) ] 
    w = $op(u,v)
    test_index_cell_array(w,o)

  end
end

# Eq

u = TestIterCellValue(a,l)
v = TestIterCellValue(a,l)
@test u == v
@test u ≈ v

# cellsum

D = 3
a = rand(2,3,4)
r = reshape(sum(a,dims=D),(2,3))
o = [CachedArray(r) for i in 1:l]

u = TestIterCellValue(a,l)
w = cellsum(u,dim=D)
test_iter_cell_array(w,o)

u = TestIndexCellValue(a,l)
w = cellsum(u,dim=D)
test_index_cell_array(w,o)

a = rand(4)
r = sum(a)
o = [r for i in 1:l]

u = TestIterCellValue(a,l)
w = cellsum(u,dim=1)
@test isa(w,CellNumber)
test_iter_cell_value(w,o)

# cellmean

a = rand(2,3,4)
r = sum(a) / length(a)
o = [r for i in 1:l]

u = TestIterCellValue(a,l)
w = cellmean(u)
@test isa(w,CellNumber)
test_iter_cell_value(w,o)

# CellMaps Operations

a = VectorValue(10,10)
b = VectorValue(15,20)
p1 = VectorValue(1,1)
p2 = VectorValue(2,2)
p3 = VectorValue(3,3)
p = [p1,p2,p3]
m = MockMap(a)
r = evaluate(m,p)

for op in (:+,:-)
  @eval begin

    cm = TestIterCellValue(m,l)
    cp = TestIterCellValue(p,l)
    rm = [ CachedArray(broadcast($op,r)) for i in 1:l]
    cm2 = $op(cm)
    test_iter_cell_map(cm2,cp,rm)

    cm = TestIndexCellValue(m,l)
    cp = TestIndexCellValue(p,l)
    rm = [ CachedArray(broadcast($op,r)) for i in 1:l]
    cm2 = $op(cm)
    ca = evaluate(cm2,cp)
    test_index_cell_array(ca,rm)

  end
end

for op in (:inner,:outer)
  @eval begin

    cm = TestIterCellValue(m,l)
    cp = TestIterCellValue(p,l)
    rm = [ CachedArray(broadcast($op,r,r)) for i in 1:l]
    cm2 = $op(cm,cm)
    test_iter_cell_map(cm2,cp,rm)

    cm = TestIndexCellValue(m,l)
    cp = TestIndexCellValue(p,l)
    rm = [ CachedArray(broadcast($op,r,r)) for i in 1:l]
    cm2 = $op(cm,cm)
    ca = evaluate(cm2,cp)
    test_index_cell_array(ca,rm)

  end
end

end # module
