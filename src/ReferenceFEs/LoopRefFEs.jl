"""
    struct Loop  <: ReferenceFEName
"""
struct Loop <: ReferenceFEName end

"""
    const loop = Loop()

Singleton of the [`Loop`](@ref) reference FE name.
"""
const loop = Loop()

"""
    struct LoopConformity <: Conformity
"""
struct LoopConformity <: Conformity end
valid_conformity_symbols(::LoopConformity) = (:L2, :Loop)

"""
    LoopRefFE(::Type{V}, p::Polytope)

Reference FE for Loop splines, available on triangles.
"""
function LoopRefFE(::Type{V}, p::Polytope) where V
  D = num_dims(p)

  (is_simplex(p) && D == 2) || @unreachable "Loop reffe only defined on triangles"

  vertices = (p===TRI) ? nothing : get_vertex_coordinates(p)
  shapefuns = _box_splines_222(V, vertices)

  dofs, face_own_dofs = _loop_dof_basis(V, vertices)
  @show face_own_dofs
  conformity = LoopConformity()
  metadata = nothing

  return GenericRefFE{Loop}(
    length(shapefuns),
    p,
    shapefuns,
    dofs,
    conformity,
    metadata,
    face_own_dofs,
    shapefuns,
  )
end

function ReferenceFE(p::Polytope, ::Loop, ::Type{V}, order) where V
  order != 4 && @unreachable "Loop reffe is order 4 only."
  LoopRefFE(T,p)
end

function _box_splines_222(::Type{V}, vertices=nothing) where V
  T = eltype(V)
  prebasis = BernsteinBasisOnSimplex(Val(2), T, 4, vertices)
  _box_splines_222(V, prebasis)
end

function _box_splines_222(::Type{V}, prebasis::BernsteinBasisOnSimplex{2,T,M,4}) where {V,T,M}
  #  Bernstein polynomial multi-indices associated with:
  #   u, v, w   (barycentric coordinates)
  #   1, 2, 3   (TRI vertex/face number)
  #  α = [4, 0, 0],  Bα =  u⁴  , Bid = 1
  #  α = [3, 1, 0],  Bα = 4u³v , Bid = 2
  #  α = [3, 0, 1],  Bα = 4u³w , Bid = 3
  #  α = [2, 2, 0],  Bα = 6u²v², Bid = 4
  #  α = [2, 1, 1],  Bα =12u²vw, Bid = 5
  #  α = [2, 0, 2],  Bα = 6u²w², Bid = 6
  #  α = [1, 3, 0],  Bα = 4uv³ , Bid = 7
  #  α = [1, 2, 1],  Bα =12uv²w, Bid = 8
  #  α = [1, 1, 2],  Bα =12uvw², Bid = 9
  #  α = [1, 0, 3],  Bα = 4uw³ , Bid = 10
  #  α = [0, 4, 0],  Bα =  v⁴  , Bid = 11
  #  α = [0, 3, 1],  Bα = 4v³w , Bid = 12
  #  α = [0, 2, 2],  Bα = 6v²w², Bid = 13
  #  α = [0, 1, 3],  Bα = 4vw³ , Bid = 14
  #  α = [0, 0, 4],  Bα =  w⁴  , Bid = 15
  #
  #  References:
  #
  # Lai, M.-J.; Schumaker, L.L. (2007);
  # Spline function on Triangulations
  # https://doi.org/10.1017/CBO9780511721588
  # Figure 12.7. page 342
  #
  # Stam, J. (1998);
  # Evaluation of loop subdivision surfaces
  # Appendix A.

  scalar_change = copy(transpose(T[
#Bid 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15  # b₂₂₂ id
     2  1  0  0  0  0  0  0  0  0  0  0  0  0  0; #      1
     2  0  1  0  0  0  0  0  0  0  0  0  0  0  0; #      2
     2  3  1  4  1  0  3  1  0  0  2  1  0  0  0; #      3
    12 12 12  8 10  8  4  6  6  4  2  3  4  3  2; #      4
     2  1  3  0  1  4  0  0  1  3  0  0  0  1  2; #      5
     0  0  0  0  0  0  1  0  0  0  2  0  0  0  0; #      6
     2  4  3  8  6  4 12 10  6  3 12 12  8  4  2; #      7
     2  3  4  4  6  8  3  6 10 12  2  4  8 12 12; #      8
     0  0  0  0  0  0  0  0  0  1  0  0  0  0  2; #      9
     0  0  0  0  0  0  0  0  0  0  2  1  0  0  0; #      10
     0  0  0  0  0  0  1  1  1  1  2  3  4  3  2; #      11
     0  0  0  0  0  0  0  0  0  0  0  0  0  1  2; #      12
  ]))./24

  L = num_indep_components(V)
  if T === V
    change = scalar_change
  else
    change = zeros(V, 15, L*12)
    @inbounds for i in 1:15
      k = 1
      for j in 1:12
        s = scalar_change[i, j]
        k = Polynomials._cartprod_set_value!(change, i, s, k)
      end
    end
  end

  linear_combination(change, prebasis)
end

function _loop_dof_basis(::Type{T}, ::Nothing) where T
  _loop_dof_basis(T, get_vertex_coordinates(TRI))
end

"""
    _loop_dof_basis(::Type{T}, vertices)

Node numbering of the Loop patch, `vertices` contain the coordinates of 4,7,8.

```plain
       (1)-------(2)
      /   \\     /  \\
     /     \\   /    \\
    (3)-----(4)------(5)
   /   \\   /   \\   /   \\
  /     \\ /     \\ /     \\
(6)-----(7)-----(8)-----(9)
  \\    /  \\    /  \\    /
   \\  /    \\  /    \\  /
   (10)----(11)----(12)
```
"""
function _loop_dof_basis(::Type{T}, vertices) where T
  v4, v7, v8 = vertices
  P = eltype(vertices)
  nodes = P[         # node # barycentric coords. w.r.t. nodes (4,7,8)
    2v4 + 0v7 - 1v8, # 1    # (2, 0,-1)
    2v4 - 1v7 + 0v8, # 2    # (2,-1, 0)
    1v4 + 1v7 - 1v8, # 3    # (1, 1,-1)
    0v4 + 1v7 + 0v8, # 4    # (0, 1, 0)
    1v4 - 1v7 + 1v8, # 5    # (1,-1, 1)
    0v4 + 2v7 - 1v8, # 6    # (0, 2,-1)
    0v4 + 1v7 + 0v8, # 7    # (0, 1, 0)
    0v4 + 0v7 + 1v8, # 8    # (0, 0, 1)
    0v4 - 1v7 + 2v8, # 9    # (0,-1, 2)
   -1v4 + 2v7 + 0v8, # 10   # (1, 2, 0)
   -1v4 + 1v7 + 1v8, # 11   # (1, 1, 1)
   -1v4 + 0v7 + 2v8, # 12   # (1, 0, 2)
  ]

  dofs = LagrangianDofBasis(T, nodes; node_major=false)

  face_own_dofs = Vector{Int}[
    findall(==(4), dofs.dof_to_node),
    findall(==(7), dofs.dof_to_node),
    findall(==(8), dofs.dof_to_node),
    [], [], [], []
  ]

  dofs, face_own_dofs
end

