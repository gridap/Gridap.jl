using Gridap

model = CartesianDiscreteModel((0, 1), 10)
V = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1), conformity=:H1)
Ω = Triangulation(model)
dx = Measure(Ω, 4)

α = Ref(1.0)
g(x) = exp(-α[] * x[1])
b(v) = ∫(g * v)dx

v = get_fe_basis(V)
cell_vec_lazy = b(v)[Ω]
cell_vec_cache = Gridap.array_cache(cell_vec_lazy)
vα1 = getindex!(cell_vec_cache, cell_vec_lazy, 1)[1]

α[] = 2.0 # The memoized values in the LazyArray cache are not valid anymore
invalidate_cache!(cell_vec_cache)
vα2 = getindex!(cell_vec_cache, cell_vec_lazy, 1)[1]
@test vα1 ≠ vα2 == b(v)[Ω][1][1]
