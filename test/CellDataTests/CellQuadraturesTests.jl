module CellQuadraturesTests

using Test
using FillArrays
using Gridap.Helpers
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Integration
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Geometry

domain = (0,1,0,1,0,1)
cells = (4,4,4)
model = simplexify(CartesianDiscreteModel(domain,cells))

trian = Triangulation(model)

degree = 2
quad = CellQuadrature(trian,degree)

x = get_cell_points(quad)
@test isa(x,CellPoint)

v = GenericCellField(get_cell_shapefuns(trian),trian,ReferenceDomain())
u = GenericCellField(lazy_map(transpose,get_cell_data(v)),v.trian,v.domain_style)

m = ∇(v)⋅∇(u)
test_array(m(x),collect(m(x)))

s = integrate(m,quad)
test_array(s,collect(s))

s = ∫( ∇(v)⋅∇(u) )*quad
test_array(s,collect(s))

s = ∫(1)*quad
@test sum(s) ≈ 1

trian_N = BoundaryTriangulation(model)
quad_N = CellQuadrature(trian_N,degree)

s = ∫( ∇(v)⋅∇(u) )*quad_N
test_array(s,collect(s))

s = ∫(1)*quad_N
@test sum(s) ≈ 6

s = ∫( x->1 )*quad_N
@test sum(s) ≈ 6

cell_measure = get_cell_measure(trian)
cell_measure_N = get_cell_measure(trian_N)
@test length(cell_measure) == num_cells(model)
@test length(cell_measure_N) == num_cells(model)
@test sum(cell_measure) ≈ 1
@test sum(cell_measure_N) ≈ 6

#using Gridap.Visualization
#writevtk(trian,"trian",celldata=["a"=>cell_measure,"b"=>cell_measure_N])




#using Gridap.Fields: MockField, MockBasis, OtherMockBasis
#
#degree = 3
#quad = Quadrature(HEX,degree)
#test_quadrature(quad)
#@test sum(get_weights(quad)) ≈ 1
#
#degree = (1,2,3)
#quad = Quadrature(HEX,degree)
#test_quadrature(quad)
#@test sum(get_weights(quad)) ≈ 1
#
#degree = 3
#quad = Quadrature(TET,degree)
#test_quadrature(quad)
#@test sum(get_weights(quad)) ≈ 0.5*1/3
#
#l = 10
#quad = CellQuadrature(degree,[QUAD,],Fill(1,l))
#q = get_coordinates(quad)
#@test isa(q,Fill)
#w = get_weights(quad)
#@test isa(w,Fill)
#
#quad = CellQuadrature(degree,[QUAD,],ones(Int,l))
#q = get_coordinates(quad)
#@test isa(q,CompressedArray)
#w = get_weights(quad)
#@test isa(w,CompressedArray)
#
#ndofs = 4
#
#z = 2.0
#v = VectorValue(3.0,1.5)
#w = VectorValue(3.4,3.5)
#a = MockBasis{2}(v,ndofs)
#b = MockBasis{2}(w,ndofs)
#c = fill(1.0,ndofs)
#f = OtherMockBasis{2}(ndofs)
#g = MockField{2}(v)
#
#l = 10
#zl = [ z for  i in 1:l]
#cl = fill(c,l)
#fl = Fill(f,l)
#ϕl = lincomb(fl,cl)
#gl = fill(g,l)
#al = Fill(a,l)
#bl = fill(b,l)
#
#gf = GenericCellField(gl,ϕl,Val(true))
#af = GenericCellField(al,ϕl,Val(true),Fill((Base.OneTo(ndofs),),l),Val((:,)))
#bf = GenericCellField(bl,ϕl,Val(true),Fill((Base.OneTo(ndofs),),l),Val((:,)))
#zf = convert_to_cell_field(zl,ϕl)
#df = af*zf
#dft = trialize_cell_basis(df)
#
#@test sum(integrate(1,ϕl,quad)) ≈ 213.33333333333323
#@test sum(integrate(zl,ϕl,quad)) ≈  426.66666666666646
#
#cm = integrate(df⋅dft,ϕl,quad)
#@test cm[1] ≈ fill(560.0,4,4)
#
#cv = integrate(df⋅v,ϕl,quad)
#@test cv[1] ≈ fill(359.9999999999998,4)
#
#quad1 = CellQuadrature(degree,[QUAD,],Fill(1,4))
#quad2 = CellQuadrature(degree,[TRI,],Fill(1,6))
#quad = lazy_append(quad1,quad2)
#@test isa(get_coordinates(quad),AppendedArray)
#@test isa(get_weights(quad),AppendedArray)
#@test isa(get_array(quad),AbstractArray{<:Quadrature})
#@test isa(get_array(quad),AppendedArray)

end # module
