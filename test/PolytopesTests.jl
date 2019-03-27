using Numa, Test
# Test to check nodes on the closure of an nface of a polytope
##
using StaticArrays
orders=[2,2]
D = length(orders)
coords = [Numa.nodescoordinates(orders[i], nodestype="Equispaced") for i=1:D]
ordsp1=orders.+1; dims = tuple(ordsp1...)
coords_arr=Array{Float64,D+1}(undef,tuple(D,dims...))
tensorfillsanti!(coords_arr,coords)
coords_arr = reshape(coords_arr,D,:)

# @santiagobadia : Can I do it with reshape ?
tpcoo = Vector{Point{D}}(undef,size(coords_arr,1))
for i in 1:size(coords_arr,1)
	tpcoo[i] = Point{D}(coords_arr[i,:])
end


@generated function tensorfillsanti!(A::Array{T,N},c) where {T,N}
	quote
		@nloops $N i A begin
		t = @ntuple $N i; d = length(t)
		dim = t[1]
		nodedim = t[dim+1]
		(@nref $N A i) = c[dim][nodedim]
	end
end
end

reinterpret(coords_arr)
size(coords)

lei = [length(coords[i]) for i in 1:D]

tpcoo = Array{Point{D},2}(undef,lei...)



for (i,j) in enumerate(tpcoo)
	println(i)
	println(j)
end

using Base.Cartesian
quote
@nloops D i tpcoo begin
	t = @ntuple D i; d = length(t)
end
end

a[] =  [ Point{D}([coords[i][i] for i in 1:D]) for j in 1:length()

p = Point{D}





c[1:D][]

m =

D=3
x = rand(3, 10)
@assert size(x,1) == 3
coords_arr = reshape( reinterpret(Float64,x), (D,size(x,2)) )
reinterpret(SVector{3,Float64}, x, (size(x,2),))

v = reinterpret(SVector{3,Float64}, m, (1000000,));

reinterpret(SVector{D,Float64}, coords_arr, (9,) );

Numa.tensorfill!(coords_arr,coords)

coords_arr = reshape( reinterpret(Float64,coords), (D,prod(ordsp1)) )



polytope = Polytope([1,1,1])
nodes = NodesArray(polytope,orders)
@test length(nodes.closurenfacenodes[end-1])==12
@test nodes.closurenfacenodes[end-1][end-1]==33
##

# Similar test on the (open for dim > 0) nface
##
orders=[2,3,2]
polytope = Polytope([1,1,1])
nodes = NodesArray(polytope,orders)
@test length(nodes.nfacenodes[end-1])==2
@test nodes.nfacenodes[end-1][end-1]==18
##

# Test to check the node coordinates
##
orders=[2,3,4]
polytope = Polytope([1,1,1])
nodes = NodesArray(polytope,orders)
@test length(nodes.coordinates)==180
@test nodes.coordinates[33,:] ≈ [1.0, 1.0/3.0, 0.0]
nfacenodes = nodes.closurenfacenodes[end-1]
coords = nodes.coordinates[nfacenodes,:]
@test (prod(coords[:,1])==1)
##

# Test to check the views of n-face set for a given n
##
polytope = Polytope([1,1,1])
for j=1:length(polytope.extrusion)+1
	for i=1:length(polytope.nfaces[polytope.dimnfs[j]])
		@test (sum(polytope.nfaces[polytope.dimnfs[j]][i].extrusion)==j-1)
	end
end
##
