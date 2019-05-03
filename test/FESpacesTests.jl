##
using Numa, Test

using Numa.Quadratures
using Numa.Polytopes
using Numa.RefFEs
# using Numa.Meshes
using Numa.FESpaces
using Numa.FESpaces: ConformingFESpace

using Numa.Polytopes: PointInt
using Numa.CellValues

using Numa.Maps
using Numa.FieldValues

using Numa.Meshes

using Numa.CellValues: CellVectorByComposition

using Numa.CellMaps

using Numa.Geometry
using Numa.Geometry.Cartesian


##

D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
phi = geomap(trian)

self = trian
coords = cellcoordinates(self)
coords[1,1]
basis = cellbasis(self)
expand(basis,coords)

##
# @santiagobadia : coords is wrong, it only iterates till cell 2 for a 2x2 trian
# Weird behaviour, it is a sub-type of AbstractArray, I can get all indices OK,
# but when I iterate, it does not work... The iterate for a AbstractArray is
# in Julia...

@show size(coords)
@show typeof(coords)
@show IndexStyle(coords)
@show coords.lid_to_gid
@show coords.gid_to_val
@show coords.cv
@show length(coords)
@show isa(coords, AbstractArray)
@show coords

for a in coords
	@show a
end

@show coords

##

order=1
orders=order*ones(Int64,D)
pol_array = celltypes(trian)
extrusion = PointInt{D}(pol_array[1])
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = ConstantCellValue(basis, ncells(trian))
quad = quadrature(trian,order=2)

fun(x::Point{2}) = x[1]^2
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun

fesp = ConformingFESpace(reffe,trian,graph)
assembler = ConformingAssembler(fesp)
funh = interpolate(fun, fesp)

funh.basis
funh.coeffs


# function interpolate(fun::Function, fesp::FESpace)
this = fesp
reffe = fesp.reffe
dofb = reffe.dofbasis
trian = fesp.trian
phi = geomap(trian)
uphys = fun ∘ phi

# @santiagobadia : Why uphys size 2 ?????

celldofs = fesp.dof_eqclass
free_dofs = zeros(Float64, fesp.num_free_dofs)

length(uphys)
length(celldofs)
freedofs = 0

for l2g in celldofs
	@show l2g
end

for imap in

for (imap,l2g) in zip(uphys,celldofs)
	# @show imap
	@show l2g
	@show evaluate(dofb,imap)
	free_dofs[l2g] = evaluate(dofb,imap)
end
shb = ConstantCellValue(reffe.shfbasis, ncells(trian))
cdofs = CellVectorFromLocalToGlobal(celldofs,free_dofs)
intu = CellFieldFromExpand(shb, cdofs)
# end



##

p_arr = []
for pi in points(grid)
  global p_arr
  p_arr = [p_arr..., pi]
end
p_arr

fun.(p_arr)

funh.coeffs.gid_to_val




evaluate(funh, ref_points) == evaluate(fun, phys_points)
evaluate(gfunh, ref_points) == evaluate(fun, phys_points)
