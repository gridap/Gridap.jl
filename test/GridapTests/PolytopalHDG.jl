
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField, Gridap.Polynomials
using Gridap.CellData, Gridap.Fields, Gridap.ReferenceFEs, Gridap.Helpers

using Gridap.FESpaces: get_facet_diameter

function get_abs_normal_vector(trian)
  function normal(c)
    t = c[2] - c[1]
    n = VectorValue(-t[2],t[1])
    n = n/norm(n)
    return n
  end
  face_coords = get_cell_coordinates(trian)
  face_normals = lazy_map(constant_field,lazy_map(normal,face_coords))
  return CellData.GenericCellField(face_normals,trian,ReferenceDomain())
end

function statically_condensed_assembly(ptopo,X,Y,M,M_test,a,l)
  # Lazily assemble the global patch-systems and 
  # perform the static condensation
  assem = FESpaces.PatchAssembler(ptopo,X,Y)
  full_matvecs = assemble_matrix_and_vector(a,l,assem,X,Y)
  sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)

  # Regular assembly of the statically-assembled systems
  patch_rows = assem.strategy.array.array[2,2].patch_rows
  patch_cols = assem.strategy.array.array[2,2].patch_cols
  matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
  matdata = ([],[],[]) # dummy matdata
  vecdata = ([],[],[]) # dummy vecdata
  data = (matvecdata, matdata, vecdata)
  A, b = assemble_matrix_and_vector(SparseMatrixAssembler(M,M_test),data)
  return A, b
end

function backward_static_condensation(ptopo,X,Y,a,l,xb)
  assem = FESpaces.PatchAssembler(ptopo,X,Y)
  patch_rows_bb = assem.strategy.array.array[2,2].patch_rows
  patch_cols_bb = assem.strategy.array.array[2,2].patch_cols

  # To discuss: can patch_rows and patch_cols be different from each other for a patch?
  if patch_rows_bb != patch_cols_bb
    @notimplemented
  else
    patchwise_xb = map(patch_rows_bb) do patch_row
      patchwise_xb = xb[patch_row]
      return patchwise_xb
    end
  end

  full_matvecs = assemble_matrix_and_vector(a,l,assem,X,Y)
  patchwise_xi_vec = lazy_map(Gridap.FESpaces.BackwardStaticCondensationMap(),full_matvecs,patchwise_xb)

  patch_rows_ii = assem.strategy.array.array[1,1].patch_rows
  patch_cols_ii = assem.strategy.array.array[1,1].patch_cols
  patchwise_xi_vecdata = ([patchwise_xi_vec,],[patch_rows_ii,],[patch_cols_ii,])

  Xi = MultiFieldFESpace([X[1],X[2]])
  Yi = MultiFieldFESpace([Y[1],Y[2]])
  assem = SparseMatrixAssembler(Xi,Yi)
  xi = assemble_vector(assem, patchwise_xi_vecdata)
  wh = FEFunction(Yi,xi)
  return wh
end

u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

nc = (16,16)
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),nc))

pmodel = Gridap.Geometry.PolytopalDiscreteModel(model)
vmodel = Gridap.Geometry.voronoi(simplexify(model))
polys = get_polytopes(vmodel)

h = maximum(map(p->get_facet_diameter(p,D), polys))

# writevtk(vmodel,"polygonal_model")

D = num_cell_dims(vmodel)
Ω = Triangulation(ReferenceFE{D}, vmodel)
Γ = Triangulation(ReferenceFE{D-1}, vmodel)

Γ = Boundary(vmodel)
grid = Γ.trian.grid

ptopo = Geometry.PatchTopology(vmodel)
Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

# Reference FEs
order = 1
# reffeV = ReferenceFE(lagrangian, VectorValue{D, Float64}, order; space=:P)
# reffeQ = ReferenceFE(lagrangian, Float64, order; space=:P)


V = FESpaces.PolytopalFESpace(Ω, VectorValue{D, Float64}, order; space=:P)
Q = FESpaces.PolytopalFESpace(Ω, Float64, order; space=:P)

reffeM = ReferenceFE(lagrangian, Float64, order; space=:P) 
M_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")
M = TrialFESpace(M_test, u)

mfs = MultiField.BlockMultiFieldStyle(2,(2,1))
Y = MultiFieldFESpace([V, Q, M_test];style=mfs)
X = MultiFieldFESpace([V, Q, M];style=mfs)

τ = 1.0 # HDG stab parameter

degree = 2*(order+1)
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)
nrel = get_normal_vector(Γp)
nabs = get_abs_normal_vector(Γp)
n = (nrel⋅nabs)⋅nabs



Πn(u) = u⋅n
Π(u) = change_domain(u,Γp,DomainStyle(u))
a((qh,uh,sh),(vh,wh,lh)) = ∫( qh⋅vh - uh*(∇⋅vh) - qh⋅∇(wh) )dΩp + ∫(sh*Πn(vh))dΓp +
                           ∫((Πn(qh) + τ*(Π(uh) - sh))*(Π(wh) + lh))dΓp
l((vh,wh,hatmh)) = ∫( f*wh )*dΩp

Asc, bsc = statically_condensed_assembly(ptopo,X,Y,M,M_test,a,l)

solver = LUSolver()
ns = numerical_setup(symbolic_setup(solver,Asc),Asc)
xb = zeros(size(bsc))
solve!(xb,ns,bsc)
sh = FEFunction(M_test,xb)

xi = backward_static_condensation(ptopo,X,Y,a,l,xb)
wh, qh = xi

