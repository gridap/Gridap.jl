module SerendipityHdivHcurlFESpacesTests

using Test
using Random
using LinearAlgebra
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays: Fill

# The same mesh with the vertices of each cell relabelled by a symmetry of the
# cube, in a scrambled order: faces of every dimension are seen in non-trivial
# vertex orders, and cells of both orientations occur.
function relabel_cells(model)
  grid = get_grid(model)
  syms = ReferenceFEs._cube_symmetries(only(get_polytopes(model)))
  perms = [s[3] for s in syms[randperm(MersenneTwister(0), length(syms))]]
  cell_nodes = get_cell_node_ids(grid)
  relabelled = [collect(cell_nodes[i])[perms[mod1(i,length(perms))]] for i in 1:length(cell_nodes)]
  x = vec(collect(get_node_coordinates(grid)))
  UnstructuredDiscreteModel(UnstructuredGrid(
    x, Table(relabelled), collect(get_reffes(grid)), collect(get_cell_type(grid))))
end

# The same mesh with its interior nodes moved: multilinear, non-affine cells
function perturb(model; δ=0.08)
  grid = get_grid(model)
  x = map(vec(collect(get_node_coordinates(grid)))) do xi
    D = length(xi)
    inner = all(c -> 1e-12 < c < 1-1e-12, Tuple(xi))
    inner ? xi + δ*VectorValue(ntuple(c -> sin(7*c + 13*sum(Tuple(xi))), D)) : xi
  end
  UnstructuredDiscreteModel(UnstructuredGrid(
    x, Table(collect.(get_cell_node_ids(grid))), collect(get_reffes(grid)), collect(get_cell_type(grid))))
end

# L² norms over the skeleton of the jumps of the normal and of the tangential
# components of a random member of the FE space
function facet_jumps(model, reffe)
  V = FESpace(model, reffe)
  uh = FEFunction(V, randn(MersenneTwister(1), num_free_dofs(V)))
  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 2*get_order(reffe)+2)
  n = get_normal_vector(Λ).⁺
  ju = uh.⁺ - uh.⁻
  jn = ju⋅n
  jt = ju - jn*n
  sqrt(sum(∫(jn*jn)dΛ)), sqrt(sum(∫(jt⋅jt)dΛ))
end

function interpolation_error(model, reffe, u)
  uh = interpolate(u, FESpace(model, reffe))
  dΩ = Measure(Triangulation(model), 2*get_order(reffe)+2)
  sqrt(sum(∫((u-uh)⋅(u-uh))dΩ))
end

square = CartesianDiscreteModel((0,1,0,1), (3,3))
cube = CartesianDiscreteModel((0,1,0,1,0,1), (2,2,2))

# In 3D, the lowest order with face DoFs
for (name, div, orders2, orders3) in (
    (sdiv, true, 1:2, 1), (trimmed_sdiv, true, 0:1, 1),
    (scurl, false, 1:2, 2), (trimmed_scurl, false, 0:1, 1))
  for (model, orders) in ((square, orders2), (cube, orders3)), o in orders
    reffe = ReferenceFE(only(get_polytopes(model)), name, Float64, o)
    relabelled = relabel_cells(model)
    models = num_cell_dims(model) == 2 ? (model, relabelled) : (relabelled,)
    for m in (models..., perturb(relabelled))
      jn, jt = facet_jumps(m, reffe)
      @test (div ? jn : jt) < 1e-12
      @test (div ? jt : jn) > 1e-2
    end

    # 𝓟ₚ(Ω;Rᴰ) is reproduced, p = o for 𝓢ₒΛᵏ and for 𝓢⁻ₒ₊₁Λᵏ, however the
    # cells are numbered, and the interpolant does not depend on the numbering
    D = num_cell_dims(model)
    u(x) = VectorValue(ntuple(c -> sum(Tuple(x))^o + c*x[c]^max(o-1,0)*x[1]^min(o,1), D))
    v(x) = VectorValue(ntuple(c -> exp(x[c])*sin(sum(Tuple(x))), D))
    @test interpolation_error(model, reffe, u) < 1e-11
    @test interpolation_error(relabelled, reffe, u) < 1e-11
    @test interpolation_error(model, reffe, v) ≈ interpolation_error(relabelled, reffe, v) rtol=1e-10
  end
end

# The change of basis: trivial on aligned 2D H(curl) meshes, and a lazy ±1
# diagonal otherwise
for model in (square, cube)
  D = num_cell_dims(model)
  relabelled = relabel_cells(model)
  for name in (sdiv, scurl)
    reffe = ReferenceFE(get_polytopes(model)[1], name, Float64, 2)
    sorted_change = FESpaces.compute_cell_bases_changes(
      name, Pushforward(typeof(name)), model, Fill(reffe,num_cells(model)), nothing)
    D == 2 && name == scurl && @test isnothing(sorted_change)

    face_signs = collect(lazy_map(FESpaces.FaceDofSignMap(model,reffe), 1:num_cells(model)))
    @test all(G -> G == I, face_signs)
    face_signs = collect(lazy_map(FESpaces.FaceDofSignMap(relabelled,reffe), 1:num_cells(relabelled)))
    @test all(G -> all(s -> abs(s) == 1, G.diag), face_signs)
    @test any(G -> any(<(0), G.diag), face_signs)

    change, invt = FESpaces.compute_cell_bases_changes(
      name, Pushforward(typeof(name)), relabelled, Fill(reffe,num_cells(relabelled)), nothing)
    @test change === invt
    @test all(G -> G isa Diagonal && all(s -> abs(s) == 1, G.diag), collect(change))
  end
end

end # module
