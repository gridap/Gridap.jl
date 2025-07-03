struct BubbleMonomialBasis{T, D, et} <: AbstractVector{Monomial}
	monomials::MonomialBasis{D, T}
	co::Vector{et}
end

"""
	BubbleMonomialBasis(::Type{T}, ::Val{D}) where {T, D}

Constructs a `BubbleMonomialBasis` on a 2D or 3D simplex of type `T`.

For a 2D triangle (TRI), the bubble function is represented as 27xy(1 - x - y),
where the factor 27 ensures that the function evaluates to 1 at the barycenter.
Similarly, for the 3D case (TET), the bubble function is 256xyz(1 - x - y - z).
"""
function BubbleMonomialBasis(::Type{T}, ::Val{2}) where {T}
	terms = [
		CartesianIndex(3, 2), # x²y 
		CartesianIndex(2, 3), # xy²
		CartesianIndex(2, 2), # xy
	]
	co = eltype(T)[-27, -27, 27]
	monomials = MonomialBasis{2}(T, (2, 2), terms)
	BubbleMonomialBasis{T, 2, eltype(T)}(monomials, co)
end

function BubbleMonomialBasis(::Type{T}, ::Val{3}) where {T}
	terms = [
		CartesianIndex(3, 2, 2), # x²yz
		CartesianIndex(2, 3, 2), # xy²z
		CartesianIndex(2, 2, 3), # xyz²
		CartesianIndex(2, 2, 2), # xyz
	]
	co = eltype(T)[-256, -256, -256, 256]
	monomials = MonomialBasis{3}(T, (2, 2, 2), terms)
	BubbleMonomialBasis{T, 3, eltype(T)}(monomials, co)
end

Base.size(::BubbleMonomialBasis{T}) where {T} = (num_components(T),)
Base.getindex(::BubbleMonomialBasis, ::Integer) = Monomial()
Base.IndexStyle(::BubbleMonomialBasis) = IndexLinear()
get_order(::BubbleMonomialBasis) = 0
return_type(::BubbleMonomialBasis{T}) where {T} = T

function _lincomb!(cache::CachedArray, vals::Matrix, b::BubbleMonomialBasis{T}) where {T}
	n = num_components(T)
	(; co) = b
	fill!(cache, zero(eltype(cache)))
	@inbounds for i ∈ axes(cache, 1)
		for j ∈ eachindex(co)
			for l ∈ axes(cache, 2)
				cache[i, l] += co[j] * vals[i, n*(j-1)+l]
			end
		end
	end
	return cache.array
end

function return_cache(b::BubbleMonomialBasis{T, D}, xs::AbstractVector{<:Point{D}}) where {T, D}
	cache = return_cache(b.monomials, xs)
	return (cache, CachedMatrix(T))
end

function evaluate!(cache, b::BubbleMonomialBasis{T, D}, xs::AbstractVector{<:Point{D}}) where {T, D}
	c1, c2 = cache
	vals = evaluate!(c1, b.monomials, xs)
	setsize!(c2, (length(xs), num_components(T)))
	r = _lincomb!(c2, vals, b)
	return r
end

function return_cache(
	fg::FieldGradientArray{1, <:BubbleMonomialBasis{T, D}},
	xs::AbstractVector{<:Point{D}}) where {T, D}

	cfg = FieldGradientArray{1}(fg.fa.monomials)
	cache = return_cache(cfg, xs)
	return (cache, CachedMatrix(gradient_type(T, zero(eltype(xs)))))
end

function evaluate!(
	cache,
	fg::FieldGradientArray{1, <:BubbleMonomialBasis{T, D}},
	xs::AbstractVector{<:Point{D}}) where {T, D}

	c1, c2 = cache
	b = fg.fa
	cfg = FieldGradientArray{1}(b.monomials)
	vals = evaluate!(c1, cfg, xs)
	setsize!(c2, (length(xs), num_components(T)))
	r = _lincomb!(c2, vals, b)
	return r
end
