#"""
#    struct Loop  <: ReferenceFEName
#"""
#struct Loop <: ReferenceFEName end
#
#"""
#    const loop = Loop()
#
#Singleton of the [`Loop`](@ref) reference FE name.
#"""
#const loop = Loop()
#
#"""
#    struct LoopConformity <: Conformity
#"""
#struct LoopConformity <: Conformity end
#valid_conformity_symbols(::LoopConformity) = (:L2, :Loop)
#
#"""
#    LoopRefFE(::Type{T}, p::Polytope, order::Integer)
#
#The `order` argument has the following meaning: the divergence of the  functions in this basis
#is in the 𝓟 space of degree `order-1`. `T` is the type of scalar components.
#"""
#function LoopRefFE(::Type{T}, p::Polytope) where T
#  D = num_dims(p)
#
#  (is_simplex(p) && D == 2) || @unreachable "Loop reffe only defined on triangles"
#
#  vertices = (p===TRI) ? nothing : get_vertex_coordinates(p)
#  prebasis = BernsteinBasisOnSimplex(Val(2), T, 4, vertices)
#
#  # /!\ the prebasis has 15 polynomials, the shapefuns have 12
#  shapefuns = _box_splines_222(T, prebasis)
#  dofs = TODO
#  face_own_dofs = Vector{Int}[[4], [7], [8], [], [], [], [1,2,3,5,6,9,10,11,12]]
#  conformity = LoopConformity()
#  metadata = nothing
#
#  return GenericRefFE{Bubble}(
#    12,
#    p,
#    prebasis,
#    dofs,
#    conformity,
#    metadata,
#    face_own_dofs,
#    shapefuns,
#  )
#end
#
#function ReferenceFE(p::Polytope, ::Loop, ::Type{T}, order) where T
#  order != 4 && @unreachable "Loop reffe is order 4 only."
#  LoopRefFE(T,p)
#end

function _box_splines_222(::Type{T}, vertices=nothing) where T
  prebasis = BernsteinBasisOnSimplex(Val(2), T, 4, vertices)
  _box_splines_222(T, prebasis)
end

function _box_splines_222(::Type{T}, prebasis::BernsteinBasisOnSimplex{2,T,M,4}) where {T,M}
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


  change = transpose(T[
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
  ]./24)

  shapefuns = linear_combination(change, prebasis)
end

function array_cache_type(a)
  typeof(_array_cache!(a))
end

