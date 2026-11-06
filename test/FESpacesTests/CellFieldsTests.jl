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
using LinearAlgebra
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

# tests for hessian rules applied to sum and product of two cell-fields

@testset "Test hessian sum and product rules" begin
  domain = (0.,1.,0.,1.)
  partition = (3,3)
  model = CartesianDiscreteModel(domain,partition)

  reffe = ReferenceFE(lagrangian,Float64,2)
  V = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags="boundary")
  U = TrialFESpace(V,x -> exp(-norm(x)))

  n = num_free_dofs(U)
  Random.seed!(123)
  uh = FEFunction(U,rand(n))
  qh = FEFunction(U,rand(n))

  p = Point(0.1,0.2)

  cf1 = ∇(uh + qh)
  @test cf1(p) ≈ ∇(uh)(p) + ∇(qh)(p)

  cf2 = ∇(uh*qh)
  @test cf2(p) ≈ (qh*∇(uh))(p) + (uh*∇(qh))(p)

  cf3 = Δ(uh + qh) # equiv to tr(∇∇(⋅)), making it an OperationCellField
  @test cf3(p) ≈ Δ(uh)(p) + Δ(qh)(p)

  cf4 = ∇∇(uh + qh)
  @test cf4(p) ≈ ∇∇(uh)(p) + ∇∇(qh)(p)

  cf5 = ∇∇(uh*qh)
  @test get_array(cf5(p)) ≈ get_array(qh(p) * ∇∇(uh)(p) + uh(p) * ∇∇(qh)(p) + ∇(qh)(p) ⊗ ∇(uh)(p) + ∇(uh)(p) ⊗ ∇(qh)(p))

  cf6 = Δ(uh*qh)
  @test cf6(p) ≈ tr(get_array(qh(p) * ∇∇(uh)(p) + uh(p) * ∇∇(qh)(p) + ∇(qh)(p) ⊗ ∇(uh)(p) + ∇(uh)(p) ⊗ ∇(qh)(p)))

  # above tests for CellPoint usage
  Ω = Triangulation(model)
  dΩ = Measure(Ω,3)
  cp = get_cell_points(dΩ)

  @test prod(cf1(cp) .≈ ∇(uh)(cp) + ∇(qh)(cp))

  @test prod(cf2(cp) .≈ (qh*∇(uh))(cp) + (uh*∇(qh))(cp))

  @test prod(cf3(cp) .≈ Δ(uh)(cp) + Δ(qh)(cp))

  @test prod(cf4(cp) .≈ ∇∇(uh)(cp) + ∇∇(qh)(cp))

  @test prod(cf5(cp) .≈ (qh * ∇∇(uh))(cp) + (uh * ∇∇(qh))(cp) + (∇(qh) ⊗ ∇(uh))(cp) + (∇(uh) ⊗ ∇(qh))(cp))

  @test prod(cf6(cp) .≈ map(i -> tr.(i), cf5(cp)))
end

@testset "Test interpolation on periodic models" begin
  model = CartesianDiscreteModel((0,1,0,1),(4,4);isperiodic=(true,true))
  pts = [VectorValue(rand(2)) for i=1:30]
  cpts = compute_cell_points_from_vector_of_points(pts, Triangulation(model), PhysicalDomain())
end


# Higher-order gradients of FE functions and bases
#
# Polynomials interpolated exactly in the FE space, on Cartesian meshes with non-square cells
# (so the affine pushforward is not the identity), against their derivatives in closed form.
# A vector polynomial is a tuple of components, each a list of (exponents => coefficient).

function dmonomial(α::NTuple{D,Int}, dirs, x) where D
  β = collect(α)
  c = 1.0
  for d in dirs
    β[d] == 0 && return 0.0
    c *= β[d]
    β[d] -= 1
  end
  c * prod(x[d]^β[d] for d in 1:D)
end

polyval(P, x, dirs = ()) = sum(c * dmonomial(α, dirs, x) for (α, c) in P; init = 0.0)

# Entries of the exact N-th gradient in Gridap's storage order: derivative indices first
# (fastest), then the value components.
exact_entries(Ps, x, N, D) = [polyval(P, x, Tuple(I)) for P in Ps for I in CartesianIndices(ntuple(_ -> D, N))]

entries(v::Number) = [v]
entries(v::MultiValue) = collect(Tuple(v))

function max_error(cf, exact, dΩ)
  q = get_cell_points(dΩ)
  vals = collect(cf(q))
  xs = collect(q.cell_phys_point)
  maximum(maximum(maximum(abs, entries(v) .- exact(x)) for (v, x) in zip(vc, xc)) for (vc, xc) in zip(vals, xs))
end

@testset "Higher-order gradients" begin

  # Lagrangian Q2, 2D: every order up to the total degree, and zero above it

  model = CartesianDiscreteModel((0,2,0,1), (4,3))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 6)
  P = [(2,2)=>1.0, (2,1)=>-2.0, (1,2)=>0.5, (2,0)=>1.5, (0,1)=>1.0]
  V = FESpace(model, ReferenceFE(lagrangian, Float64, 2))
  uh = interpolate(x -> polyval(P, x), V)
  for N in 1:5
    @test max_error(gradient(uh, Val(N)), x -> exact_entries((P,), x, N, 2), dΩ) < 1e-9
  end
  q = get_cell_points(dΩ)
  @test all(collect(gradient(uh, Val(1))(q)) .== collect(∇(uh)(q)))
  @test all(collect(gradient(uh, Val(2))(q)) .== collect(∇∇(uh)(q)))
  @test all(all(iszero, v) for v in collect(gradient(uh, Val(5))(q)))

  # Lagrangian Q2, 3D

  model3 = CartesianDiscreteModel((0,1,0,2,0,1), (2,3,2))
  dΩ3 = Measure(Triangulation(model3), 4)
  P3 = [(2,2,2)=>1.0, (2,1,0)=>-1.0, (0,1,2)=>0.5]
  V3 = FESpace(model3, ReferenceFE(lagrangian, Float64, 2))
  uh3 = interpolate(x -> polyval(P3, x), V3)
  for N in 1:3
    @test max_error(gradient(uh3, Val(N)), x -> exact_entries((P3,), x, N, 3), dΩ3) < 1e-10
  end

  # Piola-mapped bases, 2D: a Q1 vector polynomial lies in RT1 and ND1

  Pv = ([(1,1)=>1.0, (1,0)=>0.5, (0,1)=>-1.0], [(1,1)=>-0.7, (0,1)=>1.2, (0,0)=>1.0])
  for (reffe, conformity) in ((raviart_thomas, :Hdiv), (nedelec, :Hcurl))
    Vv = FESpace(model, ReferenceFE(reffe, Float64, 1); conformity)
    uv = interpolate(x -> VectorValue(polyval(Pv[1], x), polyval(Pv[2], x)), Vv)
    for N in 1:3
      @test max_error(gradient(uv, Val(N)), x -> exact_entries(Pv, x, N, 2), dΩ) < 1e-9
    end
  end

  # N-th gradients of FE bases assemble, also through jumps on the skeleton

  u = get_trial_fe_basis(V)
  v = get_fe_basis(V)
  A = assemble_matrix(∫(gradient(u, Val(3)) ⊙ gradient(v, Val(3)))dΩ, V, V)
  @test A ≈ A'
  x = get_free_dof_values(uh)
  @test x' * A * x ≈ sum(∫(gradient(uh, Val(3)) ⊙ gradient(uh, Val(3)))dΩ)

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 6)
  B = assemble_matrix(∫(jump(gradient(u, Val(3))) ⊙ jump(gradient(v, Val(3))))dΛ, V, V)
  @test norm(B) > 0
  @test B ≈ B'
  @test abs(x' * B * x) < 1e-12 * norm(B) * norm(x)^2   # a global polynomial has no jumps

end

end # module
