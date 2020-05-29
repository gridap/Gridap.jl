module PrintOpTrees

using FillArrays
using Gridap
using Gridap.Arrays: CompressedArray

io = IOBuffer()

a = zeros(3)
print_op_tree(io,a)

b = Fill(4,10)
print_op_tree(io,b)

c = CompressedArray([1.0,2.0],[1,2,1,1,1,2,2,1])
print_op_tree(io,c,showid=true)

d = apply(-,collect(c))
d = apply(+,d,d)
print_op_tree(io,d)

r = reindex(d,[4,3,3,1])
print_op_tree(io,r,showid=true)

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian = Triangulation(model)
quad = CellQuadrature(trian,2)

u(x) = x[1]
cf = CellField(u, trian)

xe = get_cell_coordinates(trian)
q = get_coordinates(quad)

uq = evaluate(cf,q)

V = FESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=1,conformity=:H1)

vh = FEFunction(V,rand(num_free_dofs(V)))

vh_q = evaluate(vh,q)
∇vh_q = evaluate(∇(vh),q)

print_op_tree(io,vh_q)
print_op_tree(io,∇vh_q)

print_op_tree(io,xe)
print_op_tree(io,q)

print_op_tree(io,uq)

data = read(seekstart(io), String)

end # module
