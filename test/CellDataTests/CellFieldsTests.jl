module CellFieldsTests

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

domain = (0,1,0,1)
cells = (3,3)
model = CartesianDiscreteModel(domain,cells)

trian = Triangulation(model)
trian_N =BoundaryTriangulation(model)
trian_D =BoundaryTriangulation(model,tags="tag_8")
trian_S =SkeletonTriangulation(model)
trian_0 =Triangulation(trian_D,Int[])

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

#np = 3
#ndofs = 4
#
#p = Point(1,2)
#x = fill(p,np)
#z = 2.0
#
#v = VectorValue(3.0,1.5)
#w = VectorValue(3.4,3.5)
#a = MockBasis{2}(v,ndofs)
#b = MockBasis{2}(w,ndofs)
#c = fill(1.0,ndofs)
#f = OtherMockBasis{2}(ndofs)
#
#g = MockField{2}(v)
#
#l = 10
#xl = Fill(x,l)
#zl = [ z for  i in 1:l]
#cl = fill(c,l)
#fl = Fill(f,l)
#ϕl = lincomb(fl,cl)
#gl = fill(g,l)
#al = Fill(a,l)
#bl = fill(b,l)
#
#gf = GenericCellField(gl,ϕl,Val(true))
#gf_x = evaluate(gf,xl)
#∇gf_x = evaluate(∇(gf),xl)
#test_cell_field(gf,xl,gf_x,grad=∇gf_x)
#
#af = GenericCellField(al,ϕl,Val(true),Fill((Base.OneTo(ndofs),),l),Val((:,)))
#bf = GenericCellField(bl,ϕl,Val(true),Fill((Base.OneTo(ndofs),),l),Val((:,)))
#zf = convert_to_cell_field(zl,ϕl)
#df = af*zf
#dft = trialize_cell_basis(df)
#
#zfr = reindex(zf,[1,4,3])
#@test length(zfr) == 3
#@test length(get_array(zfr)) == 3
#@test length(get_cell_map(zfr)) == 3
#@test length(get_cell_axes(zfr)) == 3
#
## Check memoization
#df_x1 = evaluate(df,xl)
#df_x2 = evaluate(df,xl)
#@test df_x1 === df_x2
#∇gf1 = ∇(gf)
#∇gf2 = ∇(gf)
#@test ∇gf1 === ∇gf2
#@test evaluate(∇gf1,xl) === evaluate(∇gf2,xl)
#εgf1 = ε(gf)
#εgf2 = ε(gf)
#@test εgf1 === εgf2
#@test ∇×gf === ∇×gf
#@test evaluate(∇×gf,xl) === evaluate(∇×gf,xl)
#
#@test is_test(af)
#@test is_trial(dft)
#@test is_basis(af)
#@test is_basis(dft)
#mf = af⋅dft
#@test get_metasize(mf) == (:,:)
#mf_x = evaluate(mf,xl)
#@test size(mf_x[1]) == (np,ndofs,ndofs)
#
#@test get_cell_axes(mf) == Fill((Base.OneTo(ndofs),Base.OneTo(ndofs)),l)
#
#idsL = [ i*collect(1:ndofs)  for i in 1:l]
#idsR = [ 2*i*collect(1:ndofs)  for i in 1:l]
#axesL = Fill((Base.OneTo(ndofs),),l)
#axesR = Fill((Base.OneTo(ndofs),),l)
#
#idsS = merge_cell_dofs_at_skeleton(idsL,idsR,axesL,axesR)
#@test isa(idsS,VectorOfBlockArrayCoo)
#
#afS = merge_cell_fields_at_skeleton(af,2*af)
#@test isa(afS,SkeletonCellField)
#
#afL_x = evaluate(afS.plus,xl)
#afR_x = evaluate(afS.minus,xl)
#@test isa(afL_x,VectorOfBlockArrayCoo)
#@test isa(afR_x,VectorOfBlockArrayCoo)
#
#@test isa(afS*2,SkeletonCellField)
#@test isa(afS+afS,SkeletonCellField)
#
## Checks associated with trial bases
#df = bf
#dft = trialize_cell_basis(df)
#@test is_trial(dft) == true
#cell_vals = [rand(ndofs) for  i in 1:l]
#cell_field = lincomb(df,cell_vals)
#@test isa(cell_field,CellField)
#@test get_metasize(cell_field) == ()
#cell_field_x = evaluate(cell_field,xl)
#@test isa(cell_field_x[1],AbstractVector)

end # module
