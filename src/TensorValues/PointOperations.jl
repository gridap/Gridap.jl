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
# For derivatives of point valued maps
###############################################################

_as_contravariant(p::ProductTensorValue{A,VectorValue{D}}) where {A,D} =
  ProductTensorValue{A,ContraVectorValue{Con(D)}}(p.data)

# e.g. `∂ᵢ(a∘b) = ∂ⱼa(b) ∂ᵢbʲ` contracts the value index of `∂b` with the first
# partial index of `∂a`
function dot(
  a::ProductTensorValue{A,VectorValue{D}}, b::MultiValue{Sb}) where {A,D,Sb}

  D2 = first(size(b))
  contracted_a = D2  isa Cov ? _as_contravariant(a) : a
  contracted_product(Val(1), contracted_a, b)
end

dot(a::ProductTensorValue{A,VectorValue{D}}, b::VectorValue) where {A,D} =
  contracted_product(Val(1), a, b)

function dot(
  a::ProductTensorValue{A,VectorValue{D}}, b::ProductTensorValue) where {A,D}

  D2 = first(size(b))
  contracted_a = D2 isa Cov ? _as_contravariant(a) : a
  dot(contracted_a, b)
end

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
