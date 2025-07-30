module CellFieldsTests

using Test
using FillArrays
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
#using Gridap.FESpaces
using Random

domain = (0,1,0,1)
cells = (3,3)
model = CartesianDiscreteModel(domain,cells)

trian = Triangulation(model)
trian_N = BoundaryTriangulation(model)
trian_D = BoundaryTriangulation(model,tags="tag_8")
trian_S = SkeletonTriangulation(model)
trian_0 = Triangulation(model,Int32[])

trian_sv = view(trian_S,[1,5,4,2])
nv = get_normal_vector(trian_sv)
@test isa(nv,SkeletonPair)
@test isa(nv.plus,CellField)
@test isa(nv.minus,CellField)

ϕ = GenericCellField(get_cell_map(trian),trian,ReferenceDomain())
ϕinv = GenericCellField(lazy_map(inverse_map,get_cell_map(trian)),trian,PhysicalDomain())
xref = get_cell_points(trian)
xphy = CellPoint(get_array(xref),trian,PhysicalDomain())
test_array(ϕ(xref),collect1d(ϕ(xphy)),(i,j)->all(map(≈,i,j)))
test_array(ϕinv(xref),collect1d(ϕinv(xphy)),(i,j)->all(map(≈,i,j)))

x = get_cell_points(trian)
@test DomainStyle(x) == ReferenceDomain()
@test get_array(x) == get_cell_coordinates(trian)
@test get_data(x) == get_cell_ref_coordinates(trian)

px = get_physical_coordinate(trian)
test_array(px(x),collect1d(get_array(x)))

_x = change_domain(x,PhysicalDomain())
@test DomainStyle(_x) == PhysicalDomain()
@test get_array(_x) == get_cell_coordinates(trian)
@test get_data(_x) == get_cell_coordinates(trian)

_x = change_domain(x,ReferenceDomain())
@test DomainStyle(_x) == ReferenceDomain()
@test get_array(x) == get_cell_coordinates(trian)
@test get_data(x) == get_cell_ref_coordinates(trian)

ffun(x) = 2*x[1]
f = CellField(ffun,trian)
fx = f(x)
r = map(xs->ffun.(xs),get_array(x))
r = reshape(r,length(r))
test_array(fx,r,≈)

x_0 = get_cell_points(trian_0)
fx_0 = f(x_0)
test_array(fx_0,collect(fx_0))

n_S = get_normal_vector(trian_S)
x_S = get_cell_points(trian_S)

ts = get_triangulation(∇(f))
tt = get_triangulation(n_S.plus)

@test is_change_possible(ts,tt) == true

ns1 = (Operation(sqrt)((n_S⋅n_S)))
ns2 = Operation(norm)(n_S)
@test ns1.plus(x_S) == ns2.plus(x_S)
@test ns1.minus(x_S) == ns2.minus(x_S)

change_domain(n_S,ReferenceDomain(),PhysicalDomain())

nf_S = n_S⋅∇(f)

jnf_S = jump(n_S⋅∇(f))
jnfx_S = jnf_S(x_S)
test_array(jnfx_S,0*collect(jnfx_S))

h = CellField(rand(num_cells(trian_S)),trian_S)*jump(∇(f))
hx_S = h(x_S)
test_array(hx_S,collect(hx_S))

h = 3*mean(f)⋅jump(n_S⋅∇(f))
hx_S = h(x_S)
test_array(hx_S,0*collect(hx_S))

aa(f) = 4*f
test_array((aa∘f)(x),4*r,≈)

f1 = f
f2 = 2*f
b(f1,f2) = f1+f2
test_array((b∘(f1,f2))(x),3*r,≈)

f2 = ones(num_cells(trian))
test_array((b∘(f1,f2))(x),map(i->i.+1,r),≈)

v = GenericCellField(get_cell_shapefuns(trian),trian,ReferenceDomain())
vx = v(x)
test_array(vx,collect(vx))

u = GenericCellField(lazy_map(transpose,get_data(v)),v.trian,v.domain_style)
m = v*u
test_array(m(x),collect(m(x)))
m = ∇(v)⋅∇(u)
test_array(m(x),collect(m(x)))

∇vx = ∇(v)(x)
test_array(∇vx,collect(∇vx))

∇fx = ∇(f)(x)
test_array(∇fx,collect(∇fx))


k = VectorValue(1.0,2.0)
∇kfx = ((∇+k)(f))(x)
test_array(∇kfx,collect(∇kfx))

∇kvx = ((∇+k)(v))(x)
test_array(∇kvx,collect(∇kvx))

β(x) = 2*x[1]
α = CellField(x->2*x,trian)
ax = ((∇+k)(β*α))(x)
test_array(ax,collect(ax))

ν = CellField(x->2*x,trian)
ax =((∇-k)⋅ν)(x)
test_array(ax,collect(ax))

ax =((∇-k)×ν)(x)
test_array(ax,collect(ax))

ax =((∇-k)⊗ν)(x)
test_array(ax,collect(ax))

ax =(∇.*ν)(x)
test_array(ax,collect(ax))

ax =(ν.*ν)(x)
test_array(ax,collect(ax))

ax =((∇-k).*ν)(x)
test_array(ax,collect(ax))

ax =(ν⊗(∇-k))(x)
test_array(ax,collect(ax))

σ(x) = diagonal_tensor(VectorValue(1*x[1],2*x[2]))
Fields.gradient(::typeof(σ)) = x-> ThirdOrderTensorValue{2,2,2,Float64}(1,0,0,0,0,0,0,2)
ax = ((∇+k)(σ⋅α))(x)
test_array(ax,collect(ax))

h = Operation(*)(2,f)
hx = h(x)
test_array(hx,2*fx)

a = fill(2,num_cells(trian))
h = Operation(*)(a,f)
hx = h(x)
test_array(hx,2*fx)

fx = evaluate(ffun,x)
test_array(fx,r,≈)

f_N = CellField(ffun,trian_N)
x_N = get_cell_points(trian_N)
fx_N = f_N(x_N)
test_array(fx_N,collect(fx_N))

n_N = get_normal_vector(trian_N)
nx_N = n_N(x_N)
test_array(nx_N,collect(nx_N))

h = f*n_N
hx = h(x_N)
test_array(hx,collect(hx))

gfun(x) = 3*x
g = CellField(gfun,trian)

h = Operation(*)(f,g)
gx = g(x)
hx = h(x)
r = map((i,j)->broadcast(*,i,j),fx,gx)
test_array(hx,r)

h_N = Operation(*)(f_N,g)
gx_N = g(x_N)
r = map((i,j)->broadcast(*,i,j),fx_N,gx_N)
hx_N = h_N(x_N)
test_array(hx_N,r)

g_D = CellField(gfun,trian_D)

cell_h = rand(num_cells(trian))
h = CellField(cell_h,trian)
test_array(h(x),collect(h(x)))
test_array(h(x_N),collect(h(x_N)))

h_N = (2*f_N+g)⋅g
hx_N = h_N(x_N)
test_array(hx_N,collect(hx_N))

end # module
