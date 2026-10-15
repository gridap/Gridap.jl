###############################################################
# Points as contravariant vectors
###############################################################

# This file adds methods to implicitly use VectorValue/Point in place of
# ContraVectorValue in some relevant tensor operation, mainly to enable working
# with Point valued field with varianced derivatives and their chain-rule, or
# use Point/VectorValue in Differential geometry operators taking tangent
# vectors.

_as_contravariant(a::VectorValue, ::Int) = a
_as_contravariant(a::VectorValue{D}, ::Cov) where D = ContraVectorValue{D}(Tuple(a))
_as_contravariant(a::VectorValue, i::Con) = @unreachable "A point's coordinates
  are contravariant, they cannot be contracted with the contravariant index $i."

dot(a::VectorValue, b::MultiValue{Tuple{D1,D2}}) where {D1,D2} =
  contracted_product(Val(1), _as_contravariant(a,D1), b)

dot(a::MultiValue{Tuple{D1,D2}}, b::VectorValue) where {D1,D2} =
  contracted_product(Val(1), a, _as_contravariant(b,D2))

dot(a::VectorValue, b::MultiValue{Tuple{D1}}) where D1 =
  contracted_product(Val(1), _as_contravariant(a,D1), b)

dot(a::MultiValue{Tuple{D1}}, b::VectorValue) where D1 =
  contracted_product(Val(1), a, _as_contravariant(b,D1))

dot(a::VectorValue, b::VectorValue) = contracted_product(Val(1), a, b)

###############################################################
# For derivatives of point valued maps and friends
###############################################################

# The value index of the derivatives of a `Point` valued map is contravariant,
# like the coordinates it holds.
_as_contravariant(p::_PartialsTensor{O,D,VectorValue{Dv}}) where {O,D,Dv} =
  ProductTensorValue{SymNTensorValue{NTuple{O,D}},ContraVectorValue{Con(Dv)}}(p.data)

# For `∂ᵢ(a∘b) = ∂ⱼa(b) ∂ᵢbʲ` = ∂b⋅∂a(b) where `a` is a `Point` valued
# Fields, contracts the value index of `∂b` with the first partial index of
# `∂a`.
_point_partials_dot(b::_PartialsTensor{O,Dd,VectorValue{D}}, a::MultiValue) where {O,Dd,D} =
  _map_second_factor(x -> contracted_product(Val(1), _as_contravariant(x, first(size(a))), a), b)

#`∂_I(ϕ⋅b) = ∂_I(ϕ)⋅b` with `b` constant `Point`s
function dot(a::_PartialsTensor, b::VectorValue)
  _contracts_in_second_factor(typeof(a)) ||
    return contracted_product(Val(1), a, _as_contravariant(b, last(size(a))))
  _map_second_factor(x -> dot(x, b), a)
end

dot(b::_PartialsTensor{O,Dd,VectorValue{D}}, a::MultiValue) where {O,Dd,D} =
  _point_partials_dot(b, a)

dot(b::_PartialsTensor{O,Dd,VectorValue{D}}, a::VectorValue) where {O,Dd,D} =
  _point_partials_dot(b, a)

function dot(b::_PartialsTensor{O,Dd,VectorValue{D}}, a::ProductTensorValue) where {O,Dd,D}
  first(size(a)) isa Cov || return @invoke dot(b::ProductTensorValue, a::ProductTensorValue)
  dot(_as_contravariant(b), a)
end

# For `∂ᵢϕʲ = Σₐ ∂ᵢbₐ xₐʲ` where `bₐ` are scalar functions, and `xₐ` constant
# Points.
function outer(a::_PartialsTensor{O,Dd,B}, b::VectorValue{D}) where {O,Dd,B,D}
  cb = Dd isa VariancedDimension ? ContraVectorValue{D}(Tuple(b)) : b
  !iszero(num_indep_components(a)) || return contracted_product(Val(0), a, cb)
  _map_second_factor(x -> outer(x, cb), a)
end

############################
# Comparison and broadcast #
############################

(≈)(a::VectorValue{D}, b::ContraVectorValue; kwargs...) where D =
  ≈(ContraVectorValue{D}(Tuple(a)), b; kwargs...)

(≈)(a::ContraVectorValue, b::VectorValue{D}; kwargs...) where D =
  ≈(a, ContraVectorValue{D}(Tuple(b)); kwargs...)

Base.broadcasted(f, a::VectorValue{D}, b::ContraVectorValue) where D =
  Base.broadcasted(f, ContraVectorValue{D}(Tuple(a)), b)

Base.broadcasted(f, a::ContraVectorValue, b::VectorValue{D}) where D =
  Base.broadcasted(f, a, ContraVectorValue{D}(Tuple(b)))

###############################################################
# Differential geometry on un-varianced vectors
###############################################################

function interior_product(v::VectorValue{D}, ω::ExteriorFormValue) where D
  contra_v = ContraVectorValue{D}(Tuple(v))
  interior_product(contra_v, ω)
end

function flat(v::VectorValue{D}, metric_args::Vararg{_MetricArg}) where D
  contra_v = ContraVectorValue{D}(Tuple(v))
  flat(contra_v, metric_args...)
end

function flat_flat(v::VectorValue{D}) where D
  contra_v = ContraVectorValue{D}(Tuple(v))
  flat_flat(contra_v)
end

function apply_tensor(ω::ExteriorFormValue{K}, vs::Vararg{VectorValue{D},K}) where {K,D}
  K == 0 && return ω[1]
  contra_vs = map(vi -> ContraVectorValue{D}(Tuple(vi)), vs)
  apply_tensor(ω, contra_vs...)
end

function koszul(x::VectorValue{D}, ω::ExteriorFormValue) where D
  contra_x = ContraVectorValue{D}(Tuple(x))
  koszul(contra_x, ω)
end
