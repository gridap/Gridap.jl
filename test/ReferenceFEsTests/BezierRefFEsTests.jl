module BezierRefFEsTests

using Gridap

using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Fields
using Gridap.Arrays
using FillArrays

using Test

using Gridap.ReferenceFEs: _berstein_matrix
using Gridap.ReferenceFEs: berstein_basis
using Gridap.ReferenceFEs: rationalize_bernstein_basis

## Test Bernstein basis

p_filter(e,o) = sum(e) ≤ o

# 1D
p = 2
prebasis_seg = MonomialBasis{1}(Float64,p,p_filter)
C = _berstein_matrix(prebasis_seg,SEGMENT)
C12 =
[
 1  -2   1
 0   2  -2
 0   0   1
]

@test transpose(C) == C12

ϕ = berstein_basis(prebasis_seg,SEGMENT)
X = [ Point(0.0,0.0), Point(0.5,0.5), Point(1.0,0.0) ]
ψ = linear_combination(X,ϕ)
ξ = [Point(0.0),Point(0.5),Point(1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == [ Point(0.0,0.0), Point(0.5,0.25), Point(1.0,0.0) ]

p = 3
prebasis_seg = MonomialBasis{1}(Float64,p,p_filter)
C = _berstein_matrix(prebasis_seg,SEGMENT)
C13 =
[
 1  -3   3  -1
 0   3  -6   3
 0   0   3  -3
 0   0   0   1
]

@test transpose(C) == C13

ϕ = berstein_basis(prebasis_seg,SEGMENT)
X = [ Point(0.0,0.0), Point(0.25,0.5), Point(0.75,0.5), Point(1.0,0.0) ]
ψ = linear_combination(X,ϕ)
ξ = [Point(0.0),Point(0.5),Point(1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == [ Point(0.0,0.0), Point(0.5,0.375), Point(1.0,0.0) ]

# 2D

p = 2
prebasis_tri = MonomialBasis{2}(Float64,p,p_filter)
C = _berstein_matrix(prebasis_tri,TRI)
C22 =
[
 1  -2   1  -2   2   1
 0   2  -2   0  -2   0
 0   0   1   0   0   0
 0   0   0   2  -2  -2
 0   0   0   0   2   0
 0   0   0   0   0   1
]

@test transpose(C) == C22

ϕ = berstein_basis(prebasis_tri,TRI)
X = [ Point(0.0,0.0), Point(0.5,0.0), Point(1.0,0.0), Point(0.0,0.5), Point(0.5,0.5), Point(0.0,1.0) ]
ψ = linear_combination(X,ϕ)
ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(0.3,0.3),Point(1.0,0.0),Point(0.0,1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ

p = 3
prebasis_tri = MonomialBasis{2}(Float64,p,p_filter)
C = _berstein_matrix(prebasis_tri,TRI)
C23 =
[
 1  -3   3  -1  -3   6  -3   3  -3  -1
 0   3  -6   3   0  -6   6   0   3   0
 0   0   3  -3   0   0  -3   0   0   0
 0   0   0   1   0   0   0   0   0   0
 0   0   0   0   3  -6   3  -6   6   3
 0   0   0   0   0   6  -6   0  -6   0
 0   0   0   0   0   0   3   0   0   0
 0   0   0   0   0   0   0   3  -3  -3
 0   0   0   0   0   0   0   0   3   0
 0   0   0   0   0   0   0   0   0   1
]

@test transpose(C) == C23


X = [
  Point(0.0,0.0), Point(0.25,0.0), Point(0.75,0.0), Point(1.0,0.0),
  Point(0.0,0.25), Point(0.25,0.25), Point(0.75,0.25),
  Point(0.0,0.75), Point(0.25,0.75),
  Point(0.0,1.0) ]
ϕ = berstein_basis(prebasis_tri,TRI)
ψ = linear_combination(X,ϕ)
ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(0.5,0.5),Point(1.0,0.0),Point(0.0,1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ

## BezierRefFE

tri = BezierRefFE(Float64,TRI,(3,3))
nodes = get_node_coordinates(tri) * 5
ϕ = get_shapefuns(tri)
ψ = linear_combination(nodes,ϕ)
ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(0.5,0.5),Point(1.0,0.0),Point(0.0,1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ.*5

f = ConstantField(1)
q = Quadrature(TRI,3*2)
p = get_coordinates(q)
w = get_weights(q)
J = ∇(ψ)

@test integrate(f,p,w,J) ≈ (5*5)/2

tri = BezierRefFE(Float64,TRI,(3,3))
nodes = get_node_coordinates(tri)
nodes[10] = Point(0.1,0.1)
ϕ = get_shapefuns(tri)
ψ = linear_combination(nodes,ϕ)
ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(0.5,0.5),Point(1.0,0.0),Point(0.0,1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ

f = ConstantField(1)
q = Quadrature(TRI,3*2)
p = get_coordinates(q)
w = get_weights(q)
J = ∇(ψ)

@test integrate(f,p,w,J) ≈ 1/2

tri = BezierRefFE(Float64,TRI,(3,3))
nodes = get_node_coordinates(tri)
nodes[4] = Point(1/3,-0.1)
nodes[5] = Point(2/3,-0.1)
ϕ = get_shapefuns(tri)
ψ = linear_combination(nodes,ϕ)
f = ConstantField(1)
q = Quadrature(TRI,3*2)
p = get_coordinates(q)
w = get_weights(q)
J = ∇(ψ)

@test integrate(f,p,w,J) > 1/2

tet = BezierRefFE(Float64,TET,(3,3,3))
nodes = get_node_coordinates(tet) .* 5
ϕ = get_shapefuns(tet)
ψ = linear_combination(nodes,ϕ)
ξ = [Point(0.0,0.0,0.0),Point(0.0,0.0,0.5),Point(0.5,0.5,0.5),Point(1.0,0.0,0.0),Point(0.0,1.0,0.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ .* 5

f = ConstantField(1)
q = Quadrature(TET,3*2)
p = get_coordinates(q)
w = get_weights(q)
J = ∇(ψ)

@test integrate(f,p,w,J) ≈ (5*5*5)/6

quad = BezierRefFE(Float64,QUAD,(3,3))
nodes = get_node_coordinates(quad) .* 5
ϕ = get_shapefuns(quad)
ψ = linear_combination(nodes,ϕ)
ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(1.0,0.5),Point(1.0,0.0),Point(0.0,1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ .* 5

f = ConstantField(1)
q = Quadrature(QUAD,3*2)
p = get_coordinates(q)
w = get_weights(q)
J = ∇(ψ)

@test integrate(f,p,w,J) ≈ (5*5)

hex = BezierRefFE(Float64,HEX,(3,3,3))
nodes = get_node_coordinates(hex) .* 5
ϕ = get_shapefuns(hex)
ψ = linear_combination(nodes,ϕ)
ξ = [Point(0.0,0.0,0.0),Point(0.0,0.0,0.5),Point(0.5,1.0,0.5),Point(1.0,0.0,0.0),Point(0.0,1.0,0.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ .* 5

f = ConstantField(1)
q = Quadrature(HEX,3*2)
p = get_coordinates(q)
w = get_weights(q)
J = ∇(ψ)

@test integrate(f,p,w,J) ≈ (5*5*5)

tri = ReferenceFE(TRI,bezier,Float64,(3,3))
_tri = BezierRefFE(Float64,TRI,(3,3))

@test tri == _tri


tri = BezierRefFE(Float64,TRI,(3,3))
ϕ = get_shapefuns(tri)
w = ones(length(ϕ))
ϕr = rationalize_bernstein_basis(ϕ,w)
nodes = get_node_coordinates(tri)
ψ = linear_combination(nodes,ϕr)
ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(0.5,0.5),Point(1.0,0.0),Point(0.0,1.0)]
Ψ = Fill(ψ,length(ξ))
Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ

end # module
