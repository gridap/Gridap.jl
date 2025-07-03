struct Bubble <: ReferenceFEName end
const bubble = Bubble()

function BubbleRefFE(::Type{T}, p::Polytope{D}) where {T, D}
	@notimplementedif p ∉ [TRI, TET] "only TRI and TET are supported."
  x0 = mean(get_vertex_coordinates(p))
  dofs = LagrangianDofBasis(T, [x0])
	ndofs = length(dofs)
  face_dofs = [Int[] for _ in 1:num_faces(p)]
  face_dofs[end] = 1:ndofs
  prebasis = BubbleMonomialBasis(T, Val(D))
  shapefuncs = compute_shapefuns(dofs, prebasis)
  conformity = L2Conformity()
  metadata = nothing

	return GenericRefFE{Bubble}(
		ndofs,
		p,
		prebasis,
		dofs,
		conformity,
		metadata,
		face_dofs,
		shapefuncs,
	)
end

function ReferenceFE(p::Polytope, ::Bubble, ::Type{T}) where {T} 
  return BubbleRefFE(T, p)
end
