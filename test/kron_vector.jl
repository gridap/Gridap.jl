
using Gridap
using Gridap.Arrays, Gridap.Geometry, Gridap.FESpaces

model = CartesianDiscreteModel((0,1,0,1),(5,5))

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(model,reffe)

Ω = Triangulation(model)
dΩ = Measure(Ω,2)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,2)

a(u) = ∫(u)dΓ

u = get_trial_fe_basis(V)
v = get_fe_basis(V)

cu = get_array(a(u))
cv = get_array(a(v))

nu = length(cu)
nv = length(cv)
cu_ids = vcat([collect(1:nu) for i in 1:nv]...) # Fast index
cv_ids = vcat([fill(i,nu) for i in 1:nv]...)    # Slow index
dofs = get_cell_dof_ids(V)

cu_exp = CompressedArray(cu,cu_ids)
cv_exp = CompressedArray(cv,cv_ids)

uids_exp = CompressedArray(dofs,cu_ids)
vids_exp = CompressedArray(dofs,cv_ids)

δ = fill(1.0,nv*nu)

f(δ,v1,v2) = δ .* (v1 .- v2).^2
arr = lazy_map(f,δ,cv_exp,cu_exp)
matdata = ([arr],[vids_exp],[uids_exp])

assem = SparseMatrixAssembler(V,V)
assemble_matrix(assem,matdata)

