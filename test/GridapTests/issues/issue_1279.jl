module Issue1279

using Gridap
using Gridap.ReferenceFEs, Gridap.Geometry, Gridap.Arrays
using Gridap.Fields, Gridap.CellData, Gridap.FESpaces

# Build a zero-length 1D-in-2D triangulation using only Gridap types,
# replicating the internal structure of GridapEmbedded's SubFacetTriangulation.
#
# Original MWE:
#
# using Gridap, Gridap.TensorValues, Gridap.FESpaces, Gridap.Arrays
# using GridapEmbedded, GridapEmbedded.LevelSetCutters
# 
# r = 1.5 # fails (empty Γ)
# model = CartesianDiscreteModel((0,1,0,1),(3,3))
# V_φ = TestFESpace(model,ReferenceFE(lagrangian,Float64,1))
# φh = interpolate(((x,y),)->sqrt((x-0.5)^2+(y-0.5)^2)-r,V_φ)
# geo = DiscreteGeometry(φh,model)
# cutgeo = cut(model,geo)
# Γ = EmbeddedBoundary(cutgeo)
# dΓ = Measure(Γ,2)
# 
# V1 = TestFESpace(model,ReferenceFE(lagrangian,Float64,1))
# V = MultiFieldFESpace([V1,V1])
# 
# jump_u(u1,u2) = u1 - u2
# a((u1,u2),(v1,v2)) = ∫(jump_u(u1,u2)*v2)dΓ
# AffineFEOperator(a,v->0,V,V) # Fails in BroadcastingFieldOpMap
# 
# a((u1,u2),(v1,v2)) = ∫(u1*v1 + u1*v1)dΓ
# AffineFEOperator(a,v->0,V,V) # Fails in IntegrationMap

function MockSubFacetTriangulation(model)
  Dc = num_cell_dims(model)
  Dp = num_point_dims(model)

  # Empty UnstructuredGrid: 1D segments embedded in 2D space
  coords     = VectorValue{Dp,Float64}[]      # physical node coords (2D)
  rcoords    = VectorValue{Dp,Float64}[]      # reference coords in background cell (2D!)
  node_ids   = Table(Int32[], Int32[1])       # 0-row face-to-node Table
  cell_types = Int8[]
  reffe      = LagrangianRefFE(Float64, SEGMENT, 1)
  subgrid    = UnstructuredGrid(coords, node_ids, [reffe], cell_types)

  # FaceToFaceGlue: maps 1D facet ref-coords → 2D background cell ref-coords
  # This is the exact structure SubFacetTriangulation._setup_facet_ref_map builds.
  facet_to_rcoords   = lazy_map(Broadcasting(Reindex(rcoords)), node_ids)
  facet_to_shapefuns = CompressedArray([get_shapefuns(reffe)], Int8[])
  tface_to_mface_map = lazy_map(linear_combination, facet_to_rcoords, facet_to_shapefuns)
  tface_to_mface     = Int32[]
  cell_glue = FaceToFaceGlue(tface_to_mface, tface_to_mface_map, nothing)
  return GenericTriangulation(subgrid, model, (nothing, nothing, cell_glue))
end

model = CartesianDiscreteModel((0,1,0,1),(3,3))

Γ  = MockSubFacetTriangulation(model)
dΓ = Measure(Γ, 2)

V1 = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1))
V  = MultiFieldFESpace([V1, V1])
u1, u2 = get_trial_fe_basis(V)
v1, v2 = get_fe_basis(V)

# ── Failure mode 1: error during lazy_map construction (return_type) ──────────
#
# a1 fails because lazy_map(BroadcastingFieldOpMap(*), jump_block, v_block)
# calls return_type → return_value with inconsistent block sizes:
#   jump testitem: (1,1,1) from _size_zero, v testitem: (2,4) from CompressedArray
#   @check size(a,1)==size(b,1) → 1≠2 → AssertionError
#
# This error fires *during construction* of the integrand lazy array, before any
# assembly loop, so the assembler's `if length > 0` guard does not help.

a1((u1,u2),(v1,v2)) = ∫((u1 - u2)*v2)dΓ
AffineFEOperator(a1, v->0, V, V)  # FAILS: AssertionError in BroadcastingFieldOpMap

# ── Failure mode 2: error in _array_cache! (return_cache) ─────────────────────
#
# a2/a3 multifield forms do NOT fail via AffineFEOperator because the assembler
# guards `array_cache` with `if length(cellmat) > 0`.  The bug is still there:
# calling array_cache directly on the integrand array triggers the same mismatch
# inside IntegrationMap.return_cache → evaluate! → @check np==length(w)||np==0:
#   np=1 (from _size_zero), length(w)=2 (from quadrature testitem) → 1≠2≠0
#
# These are the array-level tests that expose the bug independently.

array_cache(get_array(∫(u1*v2 - u2*v2)dΓ))  # FAILS: AssertionError in IntegrationMap
array_cache(get_array(∫(u1*v1 + u1*v1)dΓ))  # FAILS: AssertionError in IntegrationMap

end # module