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

import Numa: gradient, ∇

using Numa.CellIntegration

D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
phi = geomap(trian)

order=1
orders=order*ones(Int64,D)
pol_array = celltypes(trian)
extrusion = PointInt{D}(pol_array[1])
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = ConstantCellValue(basis, ncells(trian))
quad = quadrature(trian,order=2)
gps = coordinates(quad)

fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
# fun(x::Point{2}) = x[1]^2
# gradfun(x::Point{2}) = VectorValue(2*x[1], 0.0)
# gradient(::typeof(fun)) = gradfun

fesp = ConformingFESpace(reffe,trian,graph)
assembler = ConformingAssembler(fesp)
funh = interpolate(fun, fesp)

p_arr = []
for pi in points(grid)
  global p_arr
  p_arr = [p_arr..., pi]
end
p_arr
fun.(p_arr)

@test funh.coeffs.gid_to_val == fun.(p_arr)

# Next, put negative values, etc...
D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
quad = quadrature(trian,order=4)
fun(x::Point{2}) = 1.0
ksca = integrate(fun, trian, quad)
int = sum(ksca)
@test int ≈ 1.0
# @santiagobadia : I would say that this is not the result we want...
# we want the whole integral over the whole domain ... 1.0

physbasis = attachgeomap(cellb,phi);
ab(v,u) = inner(∇(v),∇(u)) #+ inner(v,u)
am(v,u) = inner(v,u)
V = physbasis;
U = physbasis;

uphys = fun ∘ phi

evaluate(phi, gps)
evaluate(uphys, gps)
ksca = integrate(ab(uphys,uphys),trian,quad)
int = sum(ksca)
@test int ≈ 1.0
##

kvec = integrate(ab(V,uphys),trian,quad)
# @santiagobadia : Problem when showing this result, in iterate
kvec = integrate(am(V,uphys),trian,quad)
typeof(kvec)
kmat = integrate(ab(V,U),trian,quad)
# @santiagobadia : Problem in evaluation...


cellvefs = celltovefs(graph)
vefcells = veftocells(graph)

is_fixed_vef = zeros(Bool, 25)
is_fixed_vef[1:4] .= true
import Numa.FESpaces: globaldofs
gldofs, nfree, nfixed = globaldofs(reffe, cellvefs, vefcells, is_fixed_vef)
dof_eqclass = CellVectorByComposition(cellvefs, gldofs)
cellsize(dof_eqclass)
u = zeros(Float64,cellsize(dof_eqclass)...)
# v = CachedArray(u)
# setsize!(v,vsize)

# function interpolate(fun::Function, fesp::FESpace)
reffe = fesp.reffe
dofb = reffe.dofbasis
trian = fesp.trian
phi = geomap(trian)
uphys = fun ∘ phi
celldofs = fesp.dof_eqclass

celldofs = dof_eqclass
# free_dofs = zeros(Float64, fesp.num_free_dofs)
# fixed_dofs = zeros(Float64, fesp.num_fixed_dofs)
free_dofs = zeros(Float64, nfree)
fixed_dofs = zeros(Float64, nfixed)
for (imap,l2g) in zip(uphys,celldofs)
	# free_dofs[l2g] = evaluate(dofb,imap)
	b = evaluate(dofb,imap)
	for (i,gdof) in enumerate(l2g)
		if (gdof > 0)
			free_dofs[gdof] = b[i]
		else
			fixed_dofs[-gdof] = b[i]
		end
	end
end
free_dofs
fixed_dofs
shb = ConstantCellValue(reffe.shfbasis, ncells(trian))

using Numa.CellValues: CellVectorFromLocalToGlobalPosAndNeg
using Numa.CellMaps: CellFieldFromExpand
cdofs = CellVectorFromLocalToGlobalPosAndNeg(celldofs, free_dofs, fixed_dofs)
intu = CellFieldFromExpand(shb, cdofs)




# end

# function globaldofs(reffe::RefFE, cellvefs, vefcells, is_fixed_vef::AbstractVector)
# 	nfdofs=Array{Array{Int64},1}(undef,length(vefcells))
# 	c=1
# 	c_n = -1
# 	nfdofs_l = []
# 	nfdofs_g = zeros(Int, length(vefcells)+1)
# 	nfdofs_g[1] = 1
# 	for (ignf,nf) in enumerate(vefcells)
# 		owner_cell = nf[1]
# 		lid_vef = findfirst(i->i==ignf,cellvefs[owner_cell])
# 		num_nf_dofs = length(reffe.nfacedofs[lid_vef])
# 		if ( is_fixed_vef[ignf] )
# 			nfdofs_l = [nfdofs_l..., c_n:c_n-num_nf_dofs+1... ]
# 			c_n -= num_nf_dofs
# 		else
# 			nfdofs_l = [nfdofs_l..., c:c+num_nf_dofs-1... ]
# 			c += num_nf_dofs
# 		end
# 		nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
# 	end
# 	return CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g)
# end
