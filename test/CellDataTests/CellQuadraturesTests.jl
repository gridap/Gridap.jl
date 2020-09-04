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

using Gridap.Fields: MockField, MockBasis, OtherMockBasis

degree = 3
quad = Quadrature(HEX,degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

degree = (1,2,3)
quad = Quadrature(HEX,degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 1

degree = 3
quad = Quadrature(TET,degree)
test_quadrature(quad)
@test sum(get_weights(quad)) ≈ 0.5*1/3

l = 10
quad = CellQuadrature(degree,[QUAD,],Fill(1,l))
q = get_coordinates(quad)
@test isa(q,Fill)
w = get_weights(quad)
@test isa(w,Fill)

quad = CellQuadrature(degree,[QUAD,],ones(Int,l))
q = get_coordinates(quad)
@test isa(q,CompressedArray)
w = get_weights(quad)
@test isa(w,CompressedArray)

ndofs = 4

z = 2.0
v = VectorValue(3.0,1.5)
w = VectorValue(3.4,3.5)
a = MockBasis{2}(v,ndofs)
b = MockBasis{2}(w,ndofs)
c = fill(1.0,ndofs)
f = OtherMockBasis{2}(ndofs)
g = MockField{2}(v)

l = 10
zl = [ z for  i in 1:l]
cl = fill(c,l)
fl = Fill(f,l)
gl = fill(g,l)
al = Fill(a,l)
bl = fill(b,l)

ϕl = GenericCellField(lincomb(fl,cl))
gf = GenericCellField(gl)∘inverse_map(ϕl)
af = GenericCellField(al,Fill((Base.OneTo(ndofs),),l),Val((:,)))∘inverse_map(ϕl)
bf = GenericCellField(bl,Fill((Base.OneTo(ndofs),),l),Val((:,)))∘inverse_map(ϕl)
zf = convert_to_cell_field(zl,ϕl)∘inverse_map(ϕl)
df = af*zf
dft = trialize_cell_basis(df)

@test sum(integrate(1,ϕl,quad)) ≈ 213.33333333333323

@test sum(integrate(zl,ϕl,quad)) ≈  426.66666666666646

cm = integrate(df⋅dft,ϕl,quad)
@test cm[1] ≈ fill(560.0,4,4)

cv = integrate(df⋅v,ϕl,quad)
@test cv[1] ≈ fill(359.9999999999998,4)

quad1 = CellQuadrature(degree,[QUAD,],Fill(1,4))
quad2 = CellQuadrature(degree,[TRI,],Fill(1,6))
quad = lazy_append(quad1,quad2)
@test isa(get_coordinates(quad),AppendedArray)
@test isa(get_weights(quad),AppendedArray)
@test isa(get_array(quad),AbstractArray{<:Quadrature})
@test isa(get_array(quad),AppendedArray)

end # module
