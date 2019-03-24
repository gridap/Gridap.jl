
export MultivariatePolynomialBasis, evaluate!, valuetype, grad

"""
Abstract type representing a multivariate polynomial basis
with value of type T in a coordinate space of D dimensions
"""
abstract type MultivariatePolynomialBasis{D,T} end

Base.length(::MultivariatePolynomialBasis)::Int = @abstractmethod

"""
First axis of v for dofs, second for points
"""
evaluate!(::MultivariatePolynomialBasis{D,T},::Array{Point{D},1},v::Array{T,2}) where {D,T} = @abstractmethod

"""
Returns a MultivariatePolynomialBasis{TG,D} where TG
is a type whose rank is one unit grater than the one of T
"""
gradient(::MultivariatePolynomialBasis{D,T} where{D,T})::MultivariatePolynomialBasis{D,TG} = @abstractmethod

valuetype(::Type{R} where R<:MultivariatePolynomialBasis{D,T}) where {D,T} = T

valuetype(self::MultivariatePolynomialBasis) = valuetype(typeof(self))

# Concrete implementations


# TODO: This is a temporary dummy implementation that has to be deleted and
# replaced by concrete implementations that use the functionality below.
# It serves now just as an example

export ShapeFunctionsScalarQua4

struct ShapeFunctionsScalarQua4 <: MultivariatePolynomialBasis{2,Float64} end

struct GradShapeFunctionsScalarQua4 <: MultivariatePolynomialBasis{2,VectorValue{2}} end

Base.length(::ShapeFunctionsScalarQua4) = 4

function evaluate!(
  ::ShapeFunctionsScalarQua4,points::Array{Point{2},1},v::Array{Float64,2})
  for (i,point) in enumerate(points)
    xi = point[1]
    eta = point[2]
    v[1,i] = (1-xi)*(1-eta)/4.0
    v[2,i] = (1+xi)*(1-eta)/4.0
    v[3,i] = (1-xi)*(1+eta)/4.0
    v[4,i] = (1+xi)*(1+eta)/4.0
  end
end

gradient(::ShapeFunctionsScalarQua4) = GradShapeFunctionsScalarQua4()

Base.length(::GradShapeFunctionsScalarQua4) = 4

function evaluate!(
  ::GradShapeFunctionsScalarQua4,points::Array{Point{2},1},v::Array{VectorValue{2},2})
  for (i,point) in enumerate(points)
    xi = point[1]
    eta = point[2]
    v[1,i] = VectorValue{2}( (eta-1)/4.0, (xi-1)/4.0 )
    v[2,i] = VectorValue{2}( (1-eta)/4.0,-(1+xi)/4.0 )
    v[3,i] = VectorValue{2}(-(1+eta)/4.0, (1-xi)/4.0 )
    v[4,i] = VectorValue{2}( (1+eta)/4.0, (1+xi)/4.0 )
  end
end

gradient(::GradShapeFunctionsScalarQua4) = @notimplemented

# Previous functionality.
# It has to be used to implement the abstract interface above

export TensorProductPolynomialBasis
export TensorProductMonomialBasis

export UnivariatePolynomialBasis
export UnivariateMonomialBasis

export derivative, tensorproduct!, gradient, tensorproductsquare!
export gradient

"""
Abstract basis of univariate polynomials
"""
abstract type UnivariatePolynomialBasis end

"""
Univariate monomial basis of a given `order`
"""
struct UnivariateMonomialBasis <: UnivariatePolynomialBasis
    order::Int64
end

"""
Create 1-dim univariate polynomial basis of `UnivariateMonomialBasis` type
"""
function UnivariatePolynomialBasis(order::Int64)
    UnivariateMonomialBasis(order)
end

"""
Evaluate univariate monomial basis in a set of 1D points
"""
function (monomials::UnivariateMonomialBasis)(points::AbstractVector{Float64})
    dbas = monomials.order+1
    c = Array{Float64,2}(undef, dbas, length(points))
    for (j,p) ∈ enumerate(points)
        for i= 1:dbas
            c[i,j] = p^(i-1)
        end
    end
    return c
end

"""
Function to be eliminated in the future. Compute the numder-th derivative of a monomial
at a set of 1D point
"""
function derivative(monomials::UnivariateMonomialBasis, numder::Int64, x::Vector{Float64})
    c = Vector{Vector{Float64}}(undef, length(x))
    # monomials.order+1, length(x))
    for j=1:size(c,2)
        for i=1:size(c,1)
            c[i,j] = (i<=numder) ? 0.0 : prod([i-k-1 for k=0:numder-1])x[j]^(i-numder-1)
        end
    end
    return c
end

"""
Multivariate polynomial basis obtained as tensor product of univariate polynomial basis
per dimension
"""
struct TensorProductPolynomialBasis
    polynomials::Vector{UnivariatePolynomialBasis}
end

"""
Provide a `TensorProductPolynomialBasis` for a vector `order` providing the order per
dimension
"""
function TensorProductPolynomialBasis(order::Vector{Int64})
    TensorProductPolynomialBasis([UnivariatePolynomialBasis(order[i])
                                 for i=1:length(order)])
end

#"""
#TensorProductPolynomialBasis(order::Vector{Int},nodestype)
#
## Examples
#```jldoctest
#julia> a = TensorProductPolynomialBasis([2,4]; basistype="Lagrangian", nodestype="Equispaced")
#julia> b = a([0.0, 1.0])
#```
#"""

const VectorOfPoints{D} = Vector{Point{D}} where D
"""
Evaluate a `TensorProductPolynomialBasis` at an array of `Point{D}`
"""
function (a::TensorProductPolynomialBasis)(points::VectorOfPoints{D}) where {D}
    numdims = length(a.polynomials)
    @assert numdims == D "Point dim and polynomial basis dim must be identical"
    c = Vector{Array{Float64,2}}(undef,D)
    for i ∈ 1:D
        pi = [p[i] for p ∈ points]
        c[i] = a.polynomials[i](pi)
    end
    # @santiagobadia : Unfortunately, I have to do all this to get i-th coordinate
    orders = [(a.polynomials[i].order+1) for i ∈ 1:D]
    dims = tuple(orders...)
    A = Array{Float64,D}(undef,dims)
    B = Array{Float64,2}(undef,length(A),length(points))
    for j = 1:length(points)
        A.= 1
        tensorproduct!(A,c,j)
        B[:,j] = reshape(A, length(A))
    end
    return B
end

function (a::TensorProductPolynomialBasis)(numders::Int64, x::Array{Float64,1})
    numdims = length(a.polynomials)
    @assert numdims == size(x,2) "Point dim and polynomial basis dim must be identical"
    c = [[derivative(a.polynomials[i], j, x[:,i]) for j=1:numders] for i = 1:numdims]
    orders = [(a.polynomials[i].order+1) for i = 1:numdims]
    dims = tuple(orders...)
    # A=ones(Float64,dims)
    A = Array{Float64,length(c)}(undef,dims)
    A.= 1.0
    tensorproduct!(A,c)
    A = reshape(A, length(A))
    return A
end

@generated function tensorproduct!(A::Array{Float64,N},c,ip=1) where {N}
    quote
        @nloops $N i A begin
            @nexprs $N j->( (@nref $N A i) *= c[j][i_j,ip])
        end
    end
end

@generated function tensorproductsquare!(A::Array{Array{Float64,N},M},c) where {N,M}
    quote
        @nloops $N i A begin
            B = (@nref $N A i)
            @nloops $M k B begin
                @nexprs $N j ->((@nref $N B k) *= c[j][k_j,i_j])
            end
        end
    end
end

function gradient(a::TensorProductPolynomialBasis, x::Array{Array{Float64,1},1})
    spdims = length(a.polynomials)
    c = [a.polynomials[i](x[i]) for i=1:spdims]
    dc = [derivative(a.polynomials[i], 1, x[i]) for i=1:spdims]
    ordsp1 = [(a.polynomials[i].order+1) for i=1:spdims]
    gps=[length(x[i]) for i=1:spdims]
    pldims = tuple(ordsp1...)
    vals = Array{Array{Float64,spdims},spdims}(undef,(gps...))
    grad = Array{Float64,3}(undef,prod(pldims),prod(gps),spdims)
    for k=1:spdims
        for i in eachindex(vals) # linear indexing
            vals[i]=ones(Float64,pldims)
        end
        d = [ (i==k) ? dc[i] : c[i] for i=1:spdims]
        tensorproductsquare!(vals,d)
        aux = copy(vals)
        aux = reshape(aux,prod(gps))
        aux = hcat([ reshape(aux[i],length(aux[i])) for i in 1:length(aux)]...)
        grad[:,:,k] = aux
    end
    return grad
end

function gradient(a::TensorProductPolynomialBasis, x::Array{Float64,2})
    spdims = length(a.polynomials)
    c = [a.polynomials[i](x[:,i]) for i=1:spdims]
    dc = [derivative(a.polynomials[i], 1, x[:,i]) for i=1:spdims]
    orders = [(a.polynomials[i].order+1) for i=1:spdims]
    pldims = tuple(orders...)
    vals = Array{Float64,spdims}(undef,pldims)
    grad = Array{Float64,3}(undef,length(vals),size(x,1),spdims)
    for k=1:spdims
        d = [ (i==k) ? dc[i] : c[i] for i=1:spdims]
        for j = 1:size(x,1)
            vals .= 1.0
            tensorproduct!(vals,d,j)
            grad[:,j,k] = copy(reshape(vals,length(vals)))
        end
    end
    return grad
end
