struct Bubble <: ReferenceFEName end
const bubble = Bubble()

function BubbleRefFE(::Type{T}, p::Polytope{D}; type = :mini, coeffs = nothing, terms = nothing) where {T, D}
	@notimplementedif D < 1 "Bubble reference finite elements are only defined for D >= 1."

	if isnothing(coeffs) && isnothing(terms)
		if type == :mini
			terms, coeffs = _mini_bubble_terms_and_coeffs(T, p)
		else
			@notimplemented "Only :mini is supported now; however, you may also specify `terms` and `coeffs` manually."
		end
	elseif isnothing(terms) || isnothing(coeffs)
		@error "You must specify both `terms` and `coeffs` or neither."
	end

	orders = Tuple(maximum(terms) - oneunit(eltype(terms)))
	prebasis = linear_combination(coeffs, MonomialBasis{D}(T, orders, terms))
	x0 = mean(get_vertex_coordinates(p))
	dofs = LagrangianDofBasis(T, [x0])
	ndofs = length(dofs)
	face_dofs = [Int[] for _ in 1:num_faces(p)]
	face_dofs[end] = 1:ndofs
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

get_order(reffe::GenericRefFE{Bubble}) = num_vertices(reffe.polytope)

function ReferenceFE(p::Polytope, ::Bubble, ::Type{T}; kwargs...) where {T}
	return BubbleRefFE(T, p; kwargs...)
end

function _mini_bubble_terms_and_coeffs(::Type{T}, p::Polytope{D}) where {T, D}
	et = eltype(T)
	if is_simplex(p)
		# x_1x_2...x_D(1-x_1-x_2-...-x_D)
		terms = Vector{CartesianIndex{D}}(undef, D+1)
		@inbounds for j ∈ 1:D
			# x_1x_2...x_{j-1}(x_j^2)x_{j+1}...x_D
			terms[j] = CartesianIndex(ntuple(i->i==j ? 3 : 2, Val{D}()))
		end
		# x_1x_2...x_D
		terms[D+1] = CartesianIndex(tfill(2, Val{D}()))
		coeff = fill(-one(et), D+1)
		coeff[D+1] = one(et)
	elseif is_n_cube(p)
		# x_1x_2...x_D(1-x_1)(1-x_2)...(1-x_D)
		terms = Vector{CartesianIndex{D}}(undef, 2^D)
		coeff = Vector{et}(undef, 2^D)
		offset = 1
		# Loop through all terms in the binomial expansion.
		@inbounds for n ∈ 0:D
			sign = (-one(et))^n
			for idx ∈ combinations(1:D, n)
				terms[offset] = CartesianIndex(ntuple(i -> i in idx ? 3 : 2, Val{D}()))
				coeff[offset] = sign
				offset += 1
			end
		end
	else
		@notimplemented
	end

	N = num_components(T)
	coeffs = zeros(et, length(coeff) * N, N)
	@inbounds for i ∈ axes(coeffs, 2)
		coeffs[i:N:end, i] = coeff
	end

	return terms, coeffs
end


