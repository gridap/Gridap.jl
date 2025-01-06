using Gridap
using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Fields, Gridap.FESpaces, Gridap.CellData

using FillArrays

using Meshes: viz
using CairoMakie

function PolytopalFESpace(
  trian::Triangulation{Dc}, order::Integer; T = Float64, vector_type = Vector{T}
) where Dc

  function change(p)
    v = get_vertex_coordinates(p)
    xc = mean(v)
    h = 2*maximum(norm.(v .- xc))
    J = TensorValues.diagonal_tensor(VectorValue([h for i in 1:Dc]))
    AffineField(J,xc/h)
  end

  cell_prebasis = Fill(Gridap.Polynomials.MonomialBasis{Dc}(T, order), num_cells(trian))
  cell_change = map(change, get_polytopes(trian))

  cell_shapefuns = lazy_map(Broadcasting(∘), cell_prebasis, cell_change)
  cell_dof_basis = Fill(ReferenceFEs.MockDofBasis([zero(VectorValue{Dc,Float64})]), num_cells(trian))

  fe_basis, fe_dof_basis = FESpaces.compute_cell_space(
    cell_shapefuns, cell_dof_basis, PhysicalDomain(), PhysicalDomain(), trian
  )

  cell_dof_ids, nfree = Gridap.FESpaces.compute_discontinuous_cell_dofs(
    Base.OneTo(num_cells(trian)), map(length,cell_shapefuns)
  )

  ndirichlet = 0
  dirichlet_dof_tag = Int8[]
  dirichlet_cells = Int32[]
  cell_is_dirichlet = Fill(false,num_cells(trian))
  ntags = 0

  UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dof_ids,
    fe_basis,
    fe_dof_basis,
    cell_is_dirichlet,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags,
    nothing
  )
end

model = CartesianDiscreteModel((0,1,0,1),(4,4))

pmodel = Geometry.PolytopalDiscreteModel(model)
vmodel = Geometry.voronoi(simplexify(model))
polys = get_polytopes(vmodel)
# viz(vmodel;color=1:num_cells(vmodel),showpoints=true,showsegments=true)

Ω = Triangulation(vmodel)
dΩ = Measure(Ω,2)

V = PolytopalFESpace(Ω,1)

u_exact(x) = x[1] + x[2]
a(u,v) = ∫(u⋅v)dΩ
l(v) = ∫(v*u_exact)dΩ

op = AffineFEOperator(a,l,V,V)
uh = solve(op)

eh = uh - u_exact
sum(∫(eh⋅eh)dΩ)
