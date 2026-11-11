module MomentBasedQuadraturesTests

using Test
using Gridap.Helpers, Gridap.Arrays, Gridap.Fields, Gridap.TensorValues
using Gridap.Polynomials, Gridap.ReferenceFEs

function test_exactness(quad::Quadrature{D}, ref_quad, degree, filter) where D
  x, w = get_coordinates(quad), get_weights(quad)
  xr, wr = get_coordinates(ref_quad), get_weights(ref_quad)
  for e in Iterators.product(ntuple(i -> 0:degree, D)...)
    filter(e, degree) || continue
    f(x) = prod(x.data .^ e)
    @test sum(w .* f.(x)) ≈ sum(wr .* f.(xr)) atol=1e-12
  end
end

q_filter = (e, o) -> (maximum(e) <= o)
p_filter = (e, o) -> (sum(e) <= o)

# Reference polytopes

for (p, space, filter) in ((QUAD,:Q,q_filter), (HEX,:Q,q_filter), (TRI,:P,p_filter), (TET,:P,p_filter))
  D = num_dims(p)
  for degree in 1:4
    quad = MomentBasedQuadrature(p, degree; space)
    @test isa(quad, GenericQuadrature)
    @test num_points(quad) == (space == :Q ? (degree+1)^D : binomial(degree+D,D))
    test_exactness(quad, Quadrature(p, 2*degree), degree, filter)
  end
end

# Non-convex polygons: L-shape and U-shape

Lshape = Polygon(Point{2,Float64}[(0,0),(2,0),(2,1),(1,1),(1,2),(0,2)])
Ushape = Polygon(Point{2,Float64}[(0,0),(3,0),(3,2),(2,2),(2,1),(1,1),(1,2),(0,2)])

for poly in (Lshape, Ushape)
  pmin, pmax = get_bounding_box(poly)
  for degree in 1:4
    # The Q space has total degree 2*degree, so we need a richer quadrature for the moments
    k = MomentBasedQuadratureMap(2, degree; space=:Q)
    ref_quad = Quadrature(poly, 2*degree)
    quad = evaluate(k, poly, ref_quad)
    test_exactness(quad, ref_quad, degree, q_filter)
    @test sum(get_weights(quad)) ≈ get_measure(poly)

    # Nodes are placed on the bounding box (some of them outside the polytope)
    @test all(x -> all(pmin.data .- 1e-12 .<= x.data .<= pmax.data .+ 1e-12), get_coordinates(quad))
    if poly === Lshape
      @test any(x -> x[1] > 1 && x[2] > 1, get_coordinates(quad))
    end

    # P space: the default moment quadrature is already exact
    quad = MomentBasedQuadrature(poly, degree; space=:P)
    test_exactness(quad, Quadrature(poly, 2*degree), degree, p_filter)
  end
end

# Non-Kronecker path (solves with the mass matrix of the basis at the nodes)

degree = 3
k = MomentBasedQuadratureMap(2, degree)
k_nk = MomentBasedQuadratureMap(k.poly, k.nodes, k.basis, false)
ref_quad = Quadrature(Lshape, 2*degree)
quad = evaluate(k, Lshape, ref_quad)
quad_nk = evaluate(k_nk, Lshape, ref_quad)
@test get_coordinates(quad) ≈ get_coordinates(quad_nk)
@test get_weights(quad) ≈ get_weights(quad_nk)

# Cache reuse over several polytopes

polys = [Lshape, Ushape, Lshape]
quads = [Quadrature(p, 2*degree) for p in polys]
cache = return_cache(k, first(polys), first(quads))
for (p, q) in zip(polys, quads)
  quad = evaluate!(cache, k, p, q)
  @test sum(get_weights(quad)) ≈ get_measure(p)
end

@test_throws ErrorException MomentBasedQuadratureMap(2, degree; space=:S)

end # module
