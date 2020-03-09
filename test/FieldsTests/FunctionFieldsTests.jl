using Gridap
using Gridap.Inference
using Gridap.Fields
using Test

# In this file, I want to create a FunctionField.
# The idea is to create a function that takes points
# and returns a Field value (scalar, vector, tensor).
# Next, make this function cell-wise, and consider it
# as a CellField structure.

u(x) = x

p1 = Point(1,2)
p2 = Point(2,1)
p3 = Point(4,3)
p4 = Point(6,1)
x = [p1,p2,p3,p4]

# First, I want to do some research, in order to infer
# what type of thing it returns... Without more info,
# e.g., space-dim, etc, that is not possible. I would
# need to know the type of point, e.g. Point{D,Float64}.

f = FunctionField(u)
test_field(f,x,x)

Fill(f,4)
