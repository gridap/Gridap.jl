module ConformingFESpacesTests

using Test
using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields
using FillArrays

# testing compute_conforming_cell_dofs

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]
cell_reffe = expand_cell_data(reffes,get_cell_type(grid_topology))

face_labeling = get_face_labeling(model)
dirichlet_tags = ["tag_1","tag_6"]

trian = Triangulation(model)
cell_map = get_cell_map(trian)
conf = Conformity(testitem(cell_reffe))
cell_fe = CellFE(model,cell_reffe,conf)
@test get_cell_type(cell_fe) == get_cell_type(grid_topology)

cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
  CellConformity(cell_fe),grid_topology, face_labeling, dirichlet_tags)

r = [
  [-1,1,4,5,14,15,16,17,35],[1,2,5,6,18,19,17,20,36],[2,3,6,7,21,22,20,23,37],
  [4,5,8,9,15,24,25,26,38],[5,6,9,10,19,27,26,28,39],[6,7,10,11,22,29,28,30,40],
  [8,9,12,-2,24,-4,31,32,41],[9,10,-2,-3,27,-5,32,33,42],[10,11,-3,13,29,-6,33,34,43]]
test_array(cell_dofs,r)
@test nfree == 43
@test ndiri == 6
@test dirichlet_dof_tag == [1, 2, 2, 2, 2, 2]
@test dirichlet_cells == [1, 7, 8, 9]

order = 1
reffes = [LagrangianRefFE(VectorValue{2,Float64},p,order) for p in polytopes]
cell_reffe = expand_cell_data(reffes,get_cell_type(grid_topology))

dirichlet_components = [(true,true), (false,true)]

conf = Conformity(testitem(cell_reffe))
cell_fe = CellFE(model,cell_reffe,conf)

cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
  CellConformity(cell_fe),grid_topology, face_labeling, dirichlet_tags, dirichlet_components)

r = [
  [-1,1,7,9,-2,2,8,10],[1,3,9,11,2,4,10,12],[3,5,11,13,4,6,12,14],
  [7,9,15,17,8,10,16,18],[9,11,17,19,10,12,18,20],[11,13,19,21,12,14,20,22],
  [15,17,23,25,16,18,24,-3],[17,19,25,26,18,20,-3,-4],[19,21,26,27,20,22,-4,28]]

test_array(cell_dofs,r)
@test nfree==28
@test ndiri==4
@test dirichlet_dof_tag == [1, 1, 2, 2,]
@test dirichlet_cells == [1, 7, 8, 9]

order = 3
reffes = [LagrangianRefFE(VectorValue{2,Float64},p,order) for p in polytopes]
cell_reffe = expand_cell_data(reffes,get_cell_type(grid_topology))

dirichlet_components = [(true,true), (false,true)]
conf = Conformity(testitem(cell_reffe))
cell_fe = CellFE(model,cell_reffe,conf)

cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
  CellConformity(cell_fe), grid_topology, face_labeling, dirichlet_tags, dirichlet_components)

reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},3)

V = FESpace(model,reffe,dirichlet_tags=dirichlet_tags)
@test get_cell_is_dirichlet(V) === V.cell_is_dirichlet
test_single_field_fe_space(V)

matvecdata = []
matdata = []
vecdata = []
test_single_field_fe_space(V,matvecdata,matdata,vecdata,trian)

V = FESpace(model,reffe,dirichlet_tags=dirichlet_tags,dirichlet_masks=dirichlet_components)
test_single_field_fe_space(V)

V = FESpace(trian,reffe,dirichlet_tags=dirichlet_tags,dirichlet_masks=dirichlet_components)
test_single_field_fe_space(V)

######################

V = FESpace(trian,ReferenceFE(lagrangian,Float64,1))
cell_conformity = FESpaces.get_cell_conformity(V)
@test FESpaces.get_d_ctype_lface_dofs(cell_conformity, polytopes) == [
  [[[1],[2],[3],[4]]],
  [[[1, 2], [3, 4], [1, 3], [2, 4]]],
  [[[1, 2, 3, 4]]]
]

V = FESpace(trian,ReferenceFE(lagrangian,Float64,2))
cell_conformity = FESpaces.get_cell_conformity(V)
@test FESpaces.get_d_ctype_lface_dofs(cell_conformity, polytopes) == [
  [[[1], [2], [3], [4]]],
  [[[1, 2, 5], [3, 4, 6], [1, 3, 7], [2, 4, 8]]],
  [[[1, 2, 3, 4, 5, 6, 7, 8, 9]]]
]

##################

model = CartesianDiscreteModel((0,1,0,1),(3,3))
reffe = ReferenceFE(QUAD,lagrangian,Float64,1)

V = FESpace(model,reffe)
cell_conformity = FESpaces.get_cell_conformity(V)
@test isa(cell_conformity, FESpaces.CompressedCellConformity)
bmask = FESpaces.generate_dof_mask(V,get_face_labeling(model),"boundary")
@test sum(bmask) == num_free_dofs(V) - 4
bmask_rev = FESpaces.generate_dof_mask(V,get_face_labeling(model),"boundary",reverse=true)
@test sum(bmask_rev) == 4
@test all(bmask .== .!bmask_rev)

cell_lface_own_ldofs = collect(expand_cell_data(cell_conformity.ctype_lface_own_ldofs,cell_conformity.cell_ctype))
cell_d_num_dfaces = [[cell_conformity.d_ctype_num_dfaces[d+1][ctype] for d in 0:2] for ctype in cell_conformity.cell_ctype]
cell_conformity_gen = FESpaces.GenericCellConformity(cell_lface_own_ldofs, cell_d_num_dfaces)
@test num_cells(cell_conformity_gen) == 9
@test get_cell_type(cell_conformity_gen) == Base.OneTo(9)

# Test DOF scaling

function _freedof_value_absmax(
    reffe, one_function, D, L=1.e-4; simplex=false, scale_dof=false, use_global_meshsize=false
)

  partition = tfill(5-D,Val(D))
  stretching = if use_global_meshsize
    #no stretch if testing global mesh size
    _ -> x -> x
  else
    # The stretching helps ensuring that the FESpace with scaled DOFs keeps
    # conformity on arbitrary sizes of neighbooring elements
    L -> x -> x*(1 + 9norm(x)^2/(L^2))/10
  end

  domain = ntuple( i-> isodd(i) ? 0 : L, 2D)
  model = CartesianDiscreteModel(domain, partition, map=stretching(L))
  trian = simplex ? Triangulation(simplexify(model)) : Triangulation(model)
  global_meshsize = use_global_meshsize ? L/first(partition) : nothing
  fe_space = FESpace(trian, reffe; scale_dof, global_meshsize)
  one_fef = interpolate(one_function, fe_space)
  maxdof = maximum(abs.(get_free_dof_values(one_fef)))

  maxdof
end

function _test_FESpace_dof_scaling(reffe,one_function; dims=2:3, n_cube=true, simplex=true, tol=2)
  for D in dims
    if n_cube
      #maxdof1 = _freedof_value_absmax(reffe, one_function(D), D, 1)
      #maxdofL = _freedof_value_absmax(reffe, one_function(D), D)
      #println(D,"D ", reffe[1]," quad, ", "ratio ", round(maxdofL/maxdof1; sigdigits=true))

      maxdof1 = _freedof_value_absmax(reffe, one_function(D), D, 1; scale_dof=true)
      maxdofL = _freedof_value_absmax(reffe, one_function(D), D; scale_dof=true)
      #println(D,"D ", reffe[1]," quad, ", "ratio ", round(maxdofL/maxdof1; sigdigits=true), " (scaled)")
      @test  (maxdof1/tol < maxdofL < maxdof1*tol)

      maxdof1 = _freedof_value_absmax(reffe, one_function(D), D, 1; scale_dof=true, use_global_meshsize=true)
      maxdofL = _freedof_value_absmax(reffe, one_function(D), D; scale_dof=true, use_global_meshsize=true)
      #println(D,"D ", reffe[1]," quad, ", "ratio ", round(maxdofL/maxdof1; sigdigits=true), " (globally scaled)")
      @test  (maxdof1/tol < maxdofL < maxdof1*tol)
    end

    if simplex
      #maxdof1 = _freedof_value_absmax(reffe, one_function(D), D, 1; simplex)
      #maxdofL = _freedof_value_absmax(reffe, one_function(D), D   ; simplex)
      #println(D,"D ", reffe[1]," simplex, ", "ratio ", round(maxdofL/maxdof1; sigdigits=true))

      maxdof1 = _freedof_value_absmax(reffe, one_function(D), D, 1; simplex, scale_dof=true)
      maxdofL = _freedof_value_absmax(reffe, one_function(D), D   ; simplex, scale_dof=true)
      #println(D,"D ", reffe[1]," simplex, ", "ratio ", round(maxdofL/maxdof1; sigdigits=true), " (scaled)")
      @test  (maxdof1/tol < maxdofL < maxdof1*tol)

      maxdof1 = _freedof_value_absmax(reffe, one_function(D), D, 1; simplex, scale_dof=true, use_global_meshsize=true)
      maxdofL = _freedof_value_absmax(reffe, one_function(D), D   ; simplex, scale_dof=true, use_global_meshsize=true)
      #println(D,"D ", reffe[1]," simplex, ", "ratio ", round(maxdofL/maxdof1; sigdigits=true), " (globally scaled)")
      @test  (maxdof1/tol < maxdofL < maxdof1*tol)
    end
  end
end


# In theory, mapped k-formed scale with ~hᵏ, empirical scaling results:
# - (0/D)-form mapped: ~h⁰ (D-form same as 1-form because we don't use the broken Piola map in Gridap)
# - 1-form mapped: ~h¹
# - (D-1)-form mapped: ~hᴰ⁻¹
vec_value(D) = (D > 1 ? x->VectorValue{D}(tfill(1.0,Val(D))) : x->1.0)

reffe = ReferenceFE(nedelec, Float64, 3)
_test_FESpace_dof_scaling(reffe, vec_value)

reffe = ReferenceFE(nedelec2, Float64, 3)
_test_FESpace_dof_scaling(reffe, vec_value; n_cube=false)

reffe = ReferenceFE(raviart_thomas, Float64, 3)
_test_FESpace_dof_scaling(reffe, vec_value)
_test_FESpace_dof_scaling(reffe, vec_value, dims=4:4, n_cube=false)

reffe = ReferenceFE(bdm, Float64, 3)
_test_FESpace_dof_scaling(reffe, vec_value; n_cube=false)

# All trivial elements
#
# reffe = ReferenceFE(lagrangian, Float64, 4)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=1:4)
#
# reffe = ReferenceFE(bezier, Float64, 4)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=1:3)
#
# reffe = ReferenceFE(modalC0, Float64, 2)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=1:3, simplex=false)
# #_test_FESpace_dof_scaling(reffe, D -> (x->1+sum(x)^(3-D)), dims=1:3, simplex=false)
#
# reffe = ReferenceFE(serendipity, Float64, 4)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=1:3, simplex=false)
#
# reffe = ReferenceFE(modal_lagrangian, Float64, 4)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=1:3)
#
# reffe = ReferenceFE(modal_serendipity, Float64, 4)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=1:3, simplex=false)
#
# reffe = ReferenceFE(crouzeix_raviart, Float64, 1)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=2:2, n_cube=false)
#
# reffe = ReferenceFE(bubble, Float64, 1)
# _test_FESpace_dof_scaling(reffe, D -> (x->1.0), dims=2:2, n_cube=false)

end  # module
