
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
export TensorProductLagrangianPolynomialBasis

export PolynomialBasis
export MonomialBasis
export LagrangianPolynomialBasis

export derivative, tensorproduct!, gradient, tensorproductsquare!

export gradient

abstract type PolynomialBasis end

struct MonomialBasis <: PolynomialBasis
    order::Int64
end

struct LagrangianPolynomialBasis
    order::Int64
    nodescoordinates::Vector{Float64}
    weights::Vector{Float64}
end

#"""
#PolynomialBasis(order::Int;
#                basistype::String="Lagrangian",
#                nodestype::String="Equispaced")
#
#Create 1-dim polynomial basis of type [`Lagrangian`, `Monomial`] for a given order and a set of nodes. The nodes can be equispaced or Chebyshev nodes of second kind and take values [`Equispaced`, `Chebyshev`]
#
#See also: [`(polynomials::LagrangianPolynomialBasis)`](@ref)

# Examples
#```jldoctest
#julia> a = LagrangianPolynomialBasis(2,nodestype="Equispaced")
#```
"""
   Create 1-dim polynomial basis of type [`Lagrangian`, `Monomial`] for a given order and a set of nodes. The nodes can be equispaced or Chebyshev nodes of second kind and take values [`Equispaced`, `Chebyshev`]
"""
function PolynomialBasis(order::Int64; basistype::String="Lagrangian", nodestype::String="Equispaced")
    if basistype=="Lagrangian"
        polynomialbasis=LagrangianPolynomialBasis(order,nodestype=nodestype)
    elseif basistype=="Monomial"
        polynomialbasis=MonomialBasis(order)
    else
        @assert (0==1) "Polynomial basis type not implemented"
    end
    return polynomialbasis
end

function LagrangianPolynomialBasis(order::Int64; nodestype::String="Equispaced")
    @assert ((nodestype=="Equispaced") | (nodestype=="Chebyshev"))
    "Node type not implemented"
    ordp1 = order+1
    if nodestype == "Chebyshev"
        nodescoordinates = [cos((i-1.0)*pi/order) for i=1:ordp1]
        weights = [ ((j==1 || j==ordp1) ? wj=0.5
                    : wj=1.0 )*((-1.0)^(mod(j,2)+1)) for j=1:ordp1 ]
    elseif nodestype == "Equispaced"
        (order != 0) ? nodescoordinates = (2.0/order)*[ i-1 for i=1:ordp1].-1 : nodescoordinates = [0.0]
        weights = [binomial(order,j-1)*((-1.0)^(mod(j,2)+1)) for j=1:ordp1]
    end
    return LagrangianPolynomialBasis(order,nodescoordinates,weights)
end

#"""
#(polynomials::LagrangianPolynomialBasis)(x)
#
#Evaluate all the elements of the 1-dim Lagrangian polynomial basis at a given point x using the modified Lagrangian formula, see [The numerical stability of barycentric Lagrange interpolation, N. Higham, IMA Journal of Numerical Analysis (2004) 24, 547–556](http://www.maths.manchester.ac.uk/~higham/narep/narep440.pdf). I previously implemented a modification of the barycentric formula in [Barycentric Lagrange Interpolation, J-L Berrut and Ll. N. Trefethen, SIAM Review 46(3), pp. 501-517, 2004](https://people.maths.ox.ac.uk/trefethen/barycentric.pdf).
#
#"""
function (polynomials::LagrangianPolynomialBasis)(x)
    ordp1 = length(polynomials.nodescoordinates)
    aux = polynomials.nodescoordinates.-x
    praux = prod(aux)
    # Modified Lagrangian formula
    return [ (x ≈ polynomials.nodescoordinates[i]) ? 1.0 : praux*polynomials.weights[i]/aux[i] for i=1:ordp1 ]
end

#"""
#
#Evaluate all the elements of the 1-dim monomial basis at a given point x.
#
#See also: [`TensorProductMonomialPolynomialBasis`](@ref)
#"""
function (monomials::MonomialBasis)(x::Array{Float64,1})
    #return hcat([x.^(i-1) for i = 1:monomials.order+1]'...)
    c = Array{Float64,2}(undef,monomials.order+1,size(x,1))
    for j=1:size(c,2) # num eval points
        for i=1:size(c,1) # dim pol basis 1D
            c[i,j] = x[j]^(i-1)
        end
    end
    return c
end

#"""
#derivative(monomial::MonomialBasis,numder::Int64,x::Array{Float64,1})
#
#Compute the numder-th derivative of a monomial at point x
#"""
function derivative(monomials::MonomialBasis,numder::Int64,x::Array{Float64,1})
    #c = prod([i-j for j=0:numder-1])*x^(i-numder) for i=numder:monomials.order+1]
    c = Array{Float64,2}(undef,monomials.order+1,size(x,1))
    for j=1:size(c,2)
        for i=1:size(c,1)
            c[i,j] = (i<=numder) ? 0.0 : prod([i-k-1 for k=0:numder-1])x[j]^(i-numder-1)
        end
    end
    return c
    # return [ (j<=numder) ? 0 : prod([j-k-1 for k=0:numder-1])x^(j-numder-1) for j=1:monomial.order+1]
end


struct TensorProductPolynomialBasis
    polynomials::Array{PolynomialBasis,1}
end

#"""
#
#TensorProductLagrangianPolynomial(order::Int; basistype::String="Lagrangian", nodestype::String="Equispaced")
#
#Create a multi-dimensional tensor product polynomial basis [`Lagrangian`,`Monomial`] for a given order and a set of nodes. For a Lagrangian basis, the nodes can be equispaced or Chebyshev nodes of second kind and take values [`Equispaced`,`Chebyshev`]
#
## Examples
#```jldoctest
#julia> a = TensorProductPolynomialBasis([2,4]; basistype="Lagrangian", nodestype="Equispaced")
#julia> b = [1 2; 3 4]
#julia> mapslices(a,b,dims=[2])
#```
#"""
function TensorProductPolynomialBasis(order::Vector{Int64}; basistype="Lagrangian", nodestype="Equispaced")
    polynomials = TensorProductPolynomialBasis([PolynomialBasis(order[i],
                  basistype=basistype, nodestype=nodestype) for i=1:length(order)])
end

#"""
#TensorProductPolynomialBasis(order::Vector{Int},nodestype)
#
#Compute the n-dim tensor product polynomial basis given n 1D polynomial bases.
## Examples
#```jldoctest
#julia> a = TensorProductPolynomialBasis([2,4]; basistype="Lagrangian", nodestype="Equispaced")
#julia> b = a([0.0, 1.0])
#```
#"""
function (a::TensorProductPolynomialBasis)(x::Array{Float64,2})
    numdims = length(a.polynomials)
    @assert numdims == size(x,2) "Point dim and polynomial basis dim must be identical"
    c = [ a.polynomials[i](x[:,i]) for i = 1:numdims]
    orders = [(a.polynomials[i].order+1) for i = 1:numdims]
    dims = tuple(orders...)
    A = Array{Float64,length(c)}(undef,dims)
    B = Array{Float64,2}(undef,length(A),size(x,1))
    for j = 1:size(x,1)
        A.= 1
        tensorproduct!(A,c,j)
        B[:,j] = reshape(A, length(A))
    end
    return B
end

function (a::TensorProductPolynomialBasis)(x::Array{Array{Float64,1},1})
    numdims = length(a.polynomials)
    @assert numdims == length(x) "Point dim and polynomial basis dim must be identical"
    c = [ a.polynomials[i](x[i]) for i = 1:numdims]
    ordsp1 = [(a.polynomials[i].order+1) for i = 1:numdims]
    gps=[length(x[i]) for i=1:numdims]
    dims = tuple(ordsp1...)
    A = Array{Array{Float64,numdims},numdims}(undef,(gps...))
    for i in eachindex(A) # linear indexing
        A[i]=ones(Float64,dims)
    end
    tensorproductsquare!(A,c)
    A = reshape(A,prod(gps))
    A = hcat([ reshape(A[i],length(A[i])) for i in 1:length(A)]...)
    return A
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
