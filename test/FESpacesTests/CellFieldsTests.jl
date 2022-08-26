module CellFieldsTests2

using Test
using BlockArrays
using FillArrays
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using Random
using StaticArrays

# Test function evaluation

# Set reproducible random number seed
Random.seed!(0)
@testset "evaluating functions" for D in 1:3
    xmin = 0
    xmax = 1
    domain = repeat([xmin, xmax], D)
    ncells = 3
    partition = repeat([ncells], D)
    base_model = CartesianDiscreteModel(domain, partition)
    order = 2
    reffe = ReferenceFE(lagrangian, Float64, order)

    for M ∈ [:hypercubes, :simplices]
        model = (M == :hypercubes) ? base_model : simplexify(base_model)

        V = FESpace(model, reffe)

        coeff0 = rand(Float64)
        coeffs = rand(SVector{D,Float64})
        f(x) = coeffs ⋅ SVector(Tuple(x)) + coeff0
        # TODO: use this mechanism instead to project
        # Francesc Verdugo @fverdugo 13:11
        # a(u,v) = ∫( u*v )dΩ
        # l(v) = a(f,v)
        # Solve a fe problem with this weak form
        # See also tutorial 10, "Isotropic damage model", section "L2
        # projection", function "project"
        fh = interpolate_everywhere(f, V)
        fhcache = return_cache(fh, VectorValue(zeros(D)...))

        # Test Random points
        xs = [VectorValue(rand(D)...) for i in 1:10]
        for x in xs
            x = VectorValue(rand(D)...)
            fx = f(x)
            fhx = evaluate!(fhcache, fh, x)
            @test fhx ≈ fx
        end
        fhxs = fh(xs)
        @test fhxs ≈ f.(xs)

        nv = num_vertices(model) # Number of vertices
        nf = num_faces(model,D-1) # Number of faces
        trian = Triangulation(model)
        topo = GridTopology(model)

        pts = get_vertex_coordinates(topo) # Vertex coordinates
        face_nodes = get_face_nodes(model, D-1) # face-to-node numbering
        face_coords = lazy_map(Broadcasting(Reindex(pts)), face_nodes) # Get LazyArray of coordinates of face

        # Test a random vertex from the triangulation
        pt = pts[rand(1:nv)]
        fhpt = evaluate!(fhcache, fh, pt)
        @test fhpt .≈ f(pt)

        # Test a random point lying on the face of the polytope
        face_coord = face_coords[rand(1:nf)]
        λ = rand(length(face_coord));
        λ = (D > 1) ? λ./sum(λ) : λ
        pt = face_coord ⋅ λ # Point on the face
        fhpt = evaluate!(fhcache, fh, pt)
        @test fhpt .≈ f(pt)

        # Test with CellPoint
        # Build cell_to_fxs manually
        cache1,cache2 = fhcache
        ncells = num_cells(model)
        x_to_cell(x) = CellData._point_to_cell!(cache1, x)
        point_to_cell = map(x_to_cell, xs)
        cell_to_points, point_to_lpoint = make_inverse_table(point_to_cell, ncells)
        cell_to_xs = lazy_map(Broadcasting(Reindex(xs)), cell_to_points)
        cell_to_f = get_array(fh)
        cell_to_fxs = lazy_map(evaluate, cell_to_f, cell_to_xs)

        # Now build CellPoint with xs instead of building cell_to_xs
        cell_point_xs = compute_cell_points_from_vector_of_points(xs, trian, PhysicalDomain())
        cell_point_fxs = evaluate(fh, cell_point_xs)
        @test cell_point_fxs ≈ cell_to_fxs

    end
end

p = QUAD
D = num_dims(QUAD)
et = Float64
source_model = CartesianDiscreteModel((0,1,0,1),(2,2))

@testset "Test interpolation Lagrangian" begin
  # Lagrangian space -> Lagrangian space
  f(x) = x[1] + x[2]
  reffe = LagrangianRefFE(et, p, 1)
  V₁ = FESpace(source_model, reffe, conformity=:H1)
  fh = interpolate_everywhere(f, V₁)
  # Target Lagrangian Space
  reffe = LagrangianRefFE(et, p, 2)
  model = CartesianDiscreteModel((0,1,0,1),(4,4))
  V₂ = FESpace(model, reffe, conformity=:H1)

  ifh = Interpolable(fh)
  try
    interpolate_everywhere(fh, V₂)
  catch
    gh = interpolate_everywhere(ifh, V₂)
    pts = [VectorValue(rand(2)) for i=1:10]
    for pt in pts
      @test gh(pt) ≈ fh(pt)
    end
  end

  # VectorValued Lagrangian
  fᵥ(x) = VectorValue([x[1], x[1]+x[2]])
  reffe = ReferenceFE(lagrangian, VectorValue{2,et}, 1)
  V₁ = FESpace(source_model, reffe, conformity=:H1)
  fh = interpolate_everywhere(fᵥ, V₁)
  # Target
  reffe = ReferenceFE(lagrangian, VectorValue{2,et},  2)
  V₂ = FESpace(model, reffe, conformity=:H1)

  ifh = Interpolable(fh);
  gh = interpolate_everywhere(ifh, V₂)
  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end

  # Deformed mesh
  function map(x)
    if x[1]≈1.0 && x[2]≈1.0
      x = VectorValue(0.6,0.6)
    end
    x
  end
  model = CartesianDiscreteModel((0,1,0,1),(1,1),map=map) |> simplexify
  reffe = ReferenceFE(lagrangian,Float64,1)
  V = FESpace(model,reffe)
  f2(x) = x[1]+x[2]
  u = interpolate_everywhere(f2,V)
  x = VectorValue(0.45,0.45)
  @test_throws AssertionError u(x)
  sm=KDTreeSearch(num_nearest_vertices=2)
  ux = Interpolable(u;searchmethod=sm)(x)
  @test ux == 0.9
end

@testset "Test interpolation RT" begin
  # RT Space -> RT Space
  f(x) = VectorValue([x[1], x[2]])
  reffe = RaviartThomasRefFE(et, p, 0)
  V₁ = FESpace(source_model, reffe, conformity=:HDiv)
  fh = interpolate_everywhere(f, V₁);
  # Target RT Space
  reffe = RaviartThomasRefFE(et, p, 1)
  model = CartesianDiscreteModel((0,1,0,1),(40,40))
  V₂ = FESpace(model, reffe, conformity=:HDiv)

  ifh = Interpolable(fh)
  gh = interpolate_everywhere(ifh, V₂)
  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end
end

p = TRI
D = num_dims(TRI)
et = Float64
source_model = CartesianDiscreteModel((0,1,0,1),(2,2)) |> simplexify


@testset "Test interpolation BDM" begin
  # BDM Space -> BDM Space

  f(x) = VectorValue([x[1], x[2]])
  reffe = BDMRefFE(et, p, 1)
  V₁ = FESpace(source_model, reffe, conformity=:HDiv)
  fh = interpolate_everywhere(f, V₁);
  # Target RT Space
  reffe = RaviartThomasRefFE(et, p, 1)
  model = CartesianDiscreteModel((0,1,0,1),(40,40)) |> simplexify
  V₂ = FESpace(model, reffe, conformity=:HDiv)

  ifh = Interpolable(fh)
  gh = interpolate_everywhere(ifh, V₂)
  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end
end

end # module
