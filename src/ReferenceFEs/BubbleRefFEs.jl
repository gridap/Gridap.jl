"""
    struct Bubble <: ReferenceFEName
"""
struct Bubble <: ReferenceFEName end

"""
    const bubble = Bubble()

Singleton of the [`Bubble`](@ref) reference FE name.
"""
const bubble = Bubble()

"""
    BubbleRefFE(::Type{T}, p::Polytope{D}; type = :mini, coeffs = nothing, terms = nothing)

A Lagrangian bubble space used to enrich another element.
It contains `num_independent_components(T)` shape functions.

By default -- `type == :mini` -- this is the bubble for (ℙ1ᴰb–ℙ1) MINI element (with shape
functions of degree 3), it is available for simplices and n-cubes.

If `coeffs` and `terms` are given, the bubble shape functions are defined by the
[`linear_combination`](@ref) of the weight matrix `coeffs` with the
[`MonomialBasis`](@ref) defined by `terms`.
"""
function BubbleRefFE(::Type{T}, p::Polytope{D}; type = :mini, coeffs = nothing, terms = nothing) where {T, D}
  @notimplementedif D < 1 "Bubble reference finite elements are only defined for `D >= 1`."

  if isnothing(coeffs) && isnothing(terms)
    if type == :mini
      @notimplementedif !is_simplex(p) && !is_n_cube(p) "`:mini` element only implemented for simplices and n-cubes."
      terms, coeffs = _mini_bubble_terms_and_coeffs(T, p)
    else
      @notimplemented "Only `:mini` type is supported now; however, you may also specify `terms` and `coeffs` manually."
    end
  elseif isnothing(terms) || isnothing(coeffs)
    @unreachable "You must specify both `terms` and `coeffs` or neither."
  end

  orders = Tuple(maximum(terms) - oneunit(eltype(terms)))
  prebasis = linear_combination(coeffs, MonomialBasis(Val(D), T, orders, terms))
  x0 = mean(get_vertex_coordinates(p))
  dofs = LagrangianDofBasis(T, [x0])
  ndofs = length(dofs)

  msg =  "Wrong `terms` and/or `coeffs`, the current implementation assumes that the bubble space contains `num_indep_component(T)` shapefunctions, T=$T"
  @notimplementedif length(prebasis) != ndofs msg

  face_own_dofs = [Int[] for _ in 1:num_faces(p)]
  face_own_dofs[end] = 1:ndofs
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
    face_own_dofs,
    shapefuncs,
  )
end

function ReferenceFE(p::Polytope, ::Bubble, ::Type{T}, order::Int;
                     type = :mini, coeffs = nothing, terms = nothing) where {T}
  if  isnothing(terms)
    if type == :mini
      @notimplementedif order ≠ 1 "The MINI bubble is only implemented for order 1 (ℙ1ᴰ+b–ℙ1) element."
    end
  else
    max_orders = Tuple(maximum(Tuple(term), init=1) - 1 for term in terms)
    @check order == maximum(max_orders, init=0) "`order` is not consistent with the given monomial `terms`."
  end
  return BubbleRefFE(T, p; type, coeffs, terms)
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

  N = num_indep_components(T)
  coeffs = zeros(et, length(coeff) * N, N)
  @inbounds for i ∈ axes(coeffs, 2)
    coeffs[i:N:end, i] = coeff
  end

  return terms, coeffs
end

