module ReferenceFEInterfacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using LinearAlgebra

# Abstract API

struct MockRefFEName <: ReferenceFEName end
@test_throws ErrorException ReferenceFE(VERTEX, MockRefFEName(), 0)

struct MockRefFE{D} <: ReferenceFE{D} end
D = 0
reffe = MockRefFE{D}()
@test_throws ErrorException num_dofs(reffe)
@test_throws ErrorException get_polytope(reffe)
@test_throws ErrorException get_prebasis(reffe)
@test_throws ErrorException get_dof_basis(reffe)
@test_throws ErrorException Conformity(reffe)
@test_throws ErrorException Conformity(reffe,nothing)
@test_throws ErrorException Conformity(reffe,:L2)
@test Conformity(reffe, L2Conformity()) == L2Conformity()
@test_throws ErrorException get_face_dofs(reffe)
@test_throws ErrorException get_face_own_dofs(reffe, nothing)

@test num_dims(reffe) == D
@test num_dims(typeof(reffe)) == D
@test num_cell_dims(reffe) == D
@test num_cell_dims(typeof(reffe)) == D
@test num_point_dims(reffe) == D
@test num_point_dims(typeof(reffe)) == D

D = 2
T = Float64
order = 1
prebasis = MonomialBasis(Val(D),T,order)

polytope = QUAD
x = get_vertex_coordinates(polytope)
dofs = LagrangianDofBasis(T,x)

ndofs = 4

face_own_dofs = [[1],[2],[3],[4],Int[],Int[],Int[],Int[],Int[]]
a = [1]
b = Int[]
face_own_dofs_permutations = [[a],[a],[a],[a],[b,b],[b,b],[b,b],[b,b],fill(b,8)]
face_dofs = [[1],[2],[3],[4],[1,2],[3,4],[1,3],[2,4],[1,2,3,4]]
metadata = nothing

struct MockConformity <: Conformity end

conf = MockConformity()
reffe = GenericRefFE{typeof(conf)}(
  ndofs, polytope, prebasis, dofs, conf, metadata ,face_dofs)

function ReferenceFEs.get_face_own_dofs(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs
end

function ReferenceFEs.get_face_own_dofs_permutations(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs_permutations
end

test_reference_fe(reffe)

@test get_face_own_dofs(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4]]
@test get_face_own_dofs_permutations(reffe,L2Conformity()) == [[[]], [[]], [[]], [[]], [[]], [[]], [[]], [[]], [[1, 2, 3, 4]]]

shapefuns = get_shapefuns(reffe)

@test evaluate(shapefuns,x) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

@test get_own_dofs_permutations(reffe) == face_own_dofs_permutations[end]

# dof prebasis inversion from given shape functions
dofprebasis = dofs
shapefuns = prebasis # just to get nontrivial dofs change of basis matrix

reffe = GenericRefFE{typeof(conf)}(
  ndofs, polytope, dofprebasis, conf, metadata, face_dofs, shapefuns)

test_reference_fe(reffe)

@test get_face_own_dofs(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4]]
@test get_face_own_dofs_permutations(reffe,L2Conformity()) == [[[]], [[]], [[]], [[]], [[]], [[]], [[]], [[]], [[1, 2, 3, 4]]]
@test get_own_dofs_permutations(reffe) == face_own_dofs_permutations[end]

dofs = get_dof_basis(reffe)

@test evaluate(dofs,shapefuns) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

function test_reffe_permutations(reffe::ReferenceFE)
  conf = Conformity(reffe)
  p = get_polytope(reffe)
  φ = get_shapefuns(reffe)
  face_own_dofs = get_face_own_dofs(reffe)
  face_own_dofs_perms = get_face_own_dofs_permutations(reffe)
  face_dofs = get_face_dofs(reffe)

  # construct a geometrical map onto p
  p_reffe = LagrangianRefFE(p)
  map_basis = get_shapefuns(p_reffe)
  vertex_coords = copy(get_vertex_coordinates(p))
  Φ = linear_combination(vertex_coords, map_basis)
  Φ_inv = inverse_map(Φ)
  Jt = Broadcasting(∇)(Φ)

  # construct the pushed basis
  pushforward = Pushforward(get_name(reffe), conf)
  φ_pushed = if pushforward isa IdentityPiolaMap
    φ
  else
    Broadcasting(Operation(pushforward))(φ, Jt)
  end

  # choose some points x of p at which we will compare the shape-functions values
  r = get_order(reffe)
  x = get_nodes(get_dof_basis(LagrangianRefFE(Float64, p, r+2)))
  @assert length(x) > length(φ) "more evaluation points are needed to compare the bases"

  _flatten_rows(φx::Array{<:MultiValue}) = reinterpret(Float64,φx)
  _flatten_rows(φx::Array{<:Number}) = φx

  φx = _flatten_rows(evaluate(φ, x))
  cφ_pushed = return_cache(φ_pushed, x)
  cΦ        = return_cache(Φ,        x)
  cΦ_inv    = return_cache(Φ_inv,    x)
  #cJ = return_cache(J, x)

  n, m = length(x), length(φ)
  D = num_indep_components(value_type(φ))
  M = n*D
  @assert M == size(φx, 1)
  A = cholesky(φx'*φx)
  C = zeros(Float64, m, m)

#poly_perm = last(get_vertex_permutations(p))
  for poly_perm in get_vertex_permutations(p)
    # map from p to itself "permuted"
    permute!(vertex_coords, poly_perm) # this changes Φ and Φ_inv
    pushed_xq = evaluate!(cΦ_inv, Φ_inv, x);
    @assert pushed_xq == evaluate!(cΦ_inv, inverse_map(Φ), x)
    φpx = _flatten_rows(evaluate!(cφ_pushed, φ_pushed, pushed_xq));

    # check that φ_pushed spans the same space as φ
    B = φx'*φpx
    ldiv!(C, A, B)
    is_pushed_basis_same_space = φpx ≈ φx*C
    is_change_signed_perm = C'*C ≈ I
    @show poly_perm
    @show is_pushed_basis_same_space #@test
    @show is_change_signed_perm      #@test
    println()

    for (face, (own_dof, perms)) in enumerate(zip(face_own_dofs, face_own_dofs_perms))
      @show face, (own_dof, perms)
      # iddentify the permutation of face

      for dof in own_dof
      end
    end

    # reset Φ and Φ_inv
    invpermute!(vertex_coords, poly_perm)
  end
  println()

end

reffe = ReferenceFE(TRI, :P, 3, 0)
test_reffe_permutations(reffe)

reffe = ReferenceFE(TRI, :P, 2, 1)
test_reffe_permutations(reffe)

reffe = ReferenceFE(TRI, :P⁻, 2, 1)
test_reffe_permutations(reffe)

end # module
