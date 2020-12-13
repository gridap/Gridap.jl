module PrintOpTrees

using FillArrays
using Gridap.Arrays

io = IOBuffer()

a = zeros(3)
print_op_tree(io,a)

b = Fill(4,10)
print_op_tree(io,b)

c = CompressedArray([1.0,2.0],[1,2,1,1,1,2,2,1])
print_op_tree(io,c,showid=true)

d = lazy_map(-,collect(c))
d = lazy_map(+,d,d)
print_op_tree(io,d)

r = lazy_map(Reindex(d),[4,3,3,1])
print_op_tree(io,r,showid=true)

r = lazy_map(PosNegReindex([10,20,20],[-10,-20]),[1,2,-1,3,-2])
print_op_tree(io,r,showid=true)

t = Table([[1,2],[3,4,5],[6,7]])
print_op_tree(io,t,showid=true)

data = read(seekstart(io), String)

#domain = (0,1,0,1)
#partition = (3,3)
#model = CartesianDiscreteModel(domain,partition)
#
#trian = Triangulation(model)
#quad = CellQuadrature(trian,2)
#
#u(x) = x[1]
#cf = CellField(u, trian)
#
#xe = get_cell_coordinates(trian)
#q = get_coordinates(quad)
#
#uq = evaluate(cf,q)
#
#V = FESpace(model,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
#
#vh = FEFunction(V,rand(num_free_dofs(V)))
#
#vh_q = evaluate(vh,q)
#∇vh_q = evaluate(∇(vh),q)
#
#print_op_tree(io,vh_q)
#print_op_tree(io,∇vh_q)
#
#print_op_tree(io,xe)
#print_op_tree(io,q)
#
#print_op_tree(io,uq)
#
#data = read(seekstart(io), String)

end # module
