module ConvexifyTests

using Test
using Gridap
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: signed_area, get_reflex_faces, compute_tangent_space, convexify

# Helper for testing reflex vertices on 3D polygons (projects to 2D first)
function get_reflex_features(p::Polygon{D}) where D
  coords = collect(get_vertex_coordinates(p))
  indices = collect(1:length(coords))
  if D == 2
    return get_reflex_faces(coords, indices)
  else
    u, v = compute_tangent_space(Val(2), coords)
    origin = coords[1]
    coords2d = [Point(dot(c - origin, u), dot(c - origin, v)) for c in coords]
    return get_reflex_faces(coords2d, indices)
  end
end

# ============================================================================
# 2D Polygon Tests
# ============================================================================

# Triangle
p = Polygon([Point(0.0, 0.0), Point(1.0, 0.0), Point(0.5, 1.0)])
@test isempty(get_reflex_features(p))
result = convexify(p)
@test length(result) == 1
@test num_vertices(result[1]) == 3

# Convex quadrilateral
p = Polygon([Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0)])
@test isempty(get_reflex_features(p))
result = convexify(p)
@test length(result) == 1
@test num_vertices(result[1]) == 4

# Convex pentagon
p = Polygon([Point(cos(2π*i/5), sin(2π*i/5)) for i in 0:4])
@test isempty(get_reflex_features(p))
result = convexify(p)
@test length(result) == 1
@test num_vertices(result[1]) == 5

# Concave quadrilateral (arrow shape)
p = Polygon([Point(0.0, 0.0), Point(1.0, 0.5), Point(2.0, 0.0), Point(1.0, 2.0)])
@test length(get_reflex_features(p)) == 1
result = convexify(p)
@test length(result) == 2
for poly in result
  @test isempty(get_reflex_features(poly))
end

# Concave pentagon
p = Polygon([Point(0.0, 0.0), Point(2.0, 0.5), Point(1.0, 1.0), Point(2.0, 1.5), Point(0.0, 2.0)])
@test length(get_reflex_features(p)) == 1
result = convexify(p)
@test length(result) >= 2
for poly in result
  @test isempty(get_reflex_features(poly))
end

# L-shaped polygon
vertices = [Point(0.0, 0.0), Point(2.0, 0.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(1.0, 2.0), Point(0.0, 2.0)]
p = Polygon(vertices)
reflex = get_reflex_features(p)
@test length(reflex) == 1
@test reflex[1] == 4
result = convexify(p)
@test length(result) >= 2
for poly in result
  @test isempty(get_reflex_features(poly))
end

# U-shaped polygon
vertices = [Point(0.0, 0.0), Point(3.0, 0.0), Point(3.0, 2.0), Point(2.0, 2.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(1.0, 2.0), Point(0.0, 2.0)]
p = Polygon(vertices)
@test length(get_reflex_features(p)) == 2
result = convexify(p)
for poly in result
  @test isempty(get_reflex_features(poly))
end
original_area = abs(signed_area(vertices))
total_area = sum(poly -> abs(signed_area(collect(get_vertex_coordinates(poly)))), result)
@test total_area ≈ original_area atol=1e-10


# Star-shaped polygon
vertices = [Point(0.0, 0.0), Point(2.0, 0.0), Point(2.0, 1.0), Point(1.5, 0.5), Point(1.0, 1.0), Point(0.5, 0.5), Point(0.0, 1.0)]
p = Polygon(vertices)
@test length(get_reflex_features(p)) == 2
result = convexify(p)
@test length(result) >= 3
for poly in result
  @test isempty(get_reflex_features(poly))
end

# Hexagonal star
vertices = Point{2,Float64}[]
for i in 0:5
  push!(vertices, Point(2.0 * cos(i * π/3), 2.0 * sin(i * π/3)))
  push!(vertices, Point(1.0 * cos((i + 0.5) * π/3), 1.0 * sin((i + 0.5) * π/3)))
end
p = Polygon(vertices)
@test length(get_reflex_features(p)) == 6
result = convexify(p)
@test length(result) >= 6
for poly in result
  @test isempty(get_reflex_features(poly))
end

# Staircase polygon
vertices = [Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 1.0), Point(3.0, 1.0), Point(3.0, 2.0), Point(2.0, 2.0), Point(2.0, 3.0), Point(1.0, 3.0), Point(1.0, 4.0), Point(0.0, 4.0), Point(0.0, 3.0), Point(0.0, 1.0)]
p = Polygon(vertices)
@test length(get_reflex_features(p)) >= 3
result = convexify(p)
@test length(result) >= 4
for poly in result
  @test isempty(get_reflex_features(poly))
end
original_area = abs(signed_area(vertices))
total_area = sum(poly -> abs(signed_area(collect(get_vertex_coordinates(poly)))), result)
@test total_area ≈ original_area atol=1e-10

# Narrow corridor
vertices = [Point(0.0, 2.0), Point(0.0, 1.0), Point(1.0, 1.0), Point(1.0, 0.0), Point(1.1, 0.0), Point(1.1, 2.0)]
p = Polygon(vertices)
result = convexify(p)
for poly in result
  @test isempty(get_reflex_features(poly))
end
original_area = abs(signed_area(vertices))
total_area = sum(poly -> abs(signed_area(collect(get_vertex_coordinates(poly)))), result)
@test total_area ≈ original_area atol=1e-10

# CW-oriented polygon (should throw error)
p = Polygon([Point(0.0, 0.0), Point(0.0, 2.0), Point(1.0, 2.0), Point(1.0, 1.0), Point(2.0, 1.0), Point(2.0, 0.0)])
@test_throws AssertionError convexify(p)

# Area preservation (L-shape)
vertices = [Point(0.0, 0.0), Point(2.0, 0.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(1.0, 2.0), Point(0.0, 2.0)]
p = Polygon(vertices)
original_area = abs(signed_area(vertices))
result = convexify(p)
total_area = sum(poly -> abs(signed_area(collect(get_vertex_coordinates(poly)))), result)
@test total_area ≈ original_area atol=1e-10

# ============================================================================
# 3D Polygon Tests
# ============================================================================

# Convex polygon in z=1 plane
p = Polygon([Point(0.0, 0.0, 1.0), Point(1.0, 0.0, 1.0), Point(1.0, 1.0, 1.0), Point(0.0, 1.0, 1.0)])
result = convexify(p)
@test length(result) == 1
@test num_vertices(result[1]) == 4
for poly in result
  for v in get_vertex_coordinates(poly)
    @test v[3] ≈ 1.0
  end
end

# L-shaped polygon in z=2 plane
p = Polygon([Point(0.0, 0.0, 2.0), Point(2.0, 0.0, 2.0), Point(2.0, 1.0, 2.0), Point(1.0, 1.0, 2.0), Point(1.0, 2.0, 2.0), Point(0.0, 2.0, 2.0)])
result = convexify(p)
@test length(result) >= 2
for poly in result
  for v in get_vertex_coordinates(poly)
    @test v[3] ≈ 2.0
  end
end

# U-shaped polygon in x=0 plane
p = Polygon([Point(0.0, 0.0, 0.0), Point(0.0, 3.0, 0.0), Point(0.0, 3.0, 2.0), Point(0.0, 2.0, 2.0), Point(0.0, 2.0, 1.0), Point(0.0, 1.0, 1.0), Point(0.0, 1.0, 2.0), Point(0.0, 0.0, 2.0)])
result = convexify(p)
@test length(result) >= 2
for poly in result
  for v in get_vertex_coordinates(poly)
    @test v[1] ≈ 0.0 atol=1e-10
  end
end

# Tilted polygon (z = x + y)
p = Polygon([Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 2.0), Point(2.0, 1.0, 3.0), Point(1.0, 1.0, 2.0), Point(1.0, 2.0, 3.0), Point(0.0, 2.0, 2.0)])
result = convexify(p)
@test length(result) >= 2
for poly in result
  for v in get_vertex_coordinates(poly)
    @test v[3] ≈ v[1] + v[2] atol=1e-10
  end
end

# Vertical polygon in y=0 plane
p = Polygon([Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0), Point(2.0, 0.0, 1.0), Point(1.0, 0.0, 1.0), Point(1.0, 0.0, 2.0), Point(0.0, 0.0, 2.0)])
result = convexify(p)
@test length(result) >= 2
for poly in result
  for v in get_vertex_coordinates(poly)
    @test v[2] ≈ 0.0 atol=1e-10
  end
end

# Arbitrary plane (2x + 3y + z = 6)
on_plane(x, y) = Point(x, y, 6 - 2x - 3y)
p = Polygon([on_plane(0.0, 0.0), on_plane(1.0, 0.0), on_plane(1.0, 0.5), on_plane(0.5, 0.5), on_plane(0.5, 1.0), on_plane(0.0, 1.0)])
result = convexify(p)
@test length(result) >= 2
for poly in result
  for v in get_vertex_coordinates(poly)
    @test 2v[1] + 3v[2] + v[3] ≈ 6.0 atol=1e-10
  end
end

# Area preservation (L-shape in z=0 plane)
vertices = [Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0), Point(2.0, 1.0, 0.0), Point(1.0, 1.0, 0.0), Point(1.0, 2.0, 0.0), Point(0.0, 2.0, 0.0)]
p = Polygon(vertices)
result = convexify(p)
u, v = compute_tangent_space(Val(2), vertices)
origin = vertices[1]
total_area = sum(result) do poly
  verts3d = collect(get_vertex_coordinates(poly))
  verts2d = [Point(dot(c - origin, u), dot(c - origin, v)) for c in verts3d]
  abs(signed_area(Vector{Point{2,Float64}}(verts2d)))
end
@test total_area ≈ 3.0 atol=1e-10

end # module
