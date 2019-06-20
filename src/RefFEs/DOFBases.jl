module DOFBases

using Gridap
using Gridap.Helpers

using ..Polytopes
export DOFBasis
export numlocaldofs

import Gridap: evaluate, evaluate!

"""
Abstract DOF cell basis
"""
abstract type DOFBasis{D,T} end

"""
Evaluate the DOFs for a given polynomial basis
"""
function evaluate(this::DOFBasis{D,T},
  fields::Map{Point{D},N,T,M})::Array{Float64,M} where {D,T,N,M}
  size_arr = return_size(fields, (numlocaldofs(this),))
  aux = zeros(Float64, size_arr...)
  evaluate!(this,fields,aux)
  return aux
end

function evaluate!(this::DOFBasis{D,T},
  fields::Map{Point{D},N,T,N}, auxv::AbstractArray)::Array{Float64,N} where {D,T,N}
  @abstractmethod
end
# @santiagobadia : To replace the following ones

function evaluate!(this::DOFBasis{D,T},
  prebasis::Basis{D,T}, auxv::AbstractMatrix)::Array{Float64,2} where {D,T}
  @abstractmethod
end

function evaluate!(this::DOFBasis{D,T},
  prebasis::Field{D,T}, auxv::AbstractVector)::Vector{Float64} where {D,T}
  @abstractmethod end

"""
Lagrangian DOF basis, which consists on evaluating the polynomial basis
(prebasis) in a set of points (nodes)
"""
struct LagrangianDOFBasis{D,T} <: DOFBasis{D,T}
  nodes::Array{Point{D,Float64}}
end

"""
Evaluate the Lagrangian DOFs basis (i.e., nodal values) for a given polynomial
basis
"""
function evaluate!(this::LagrangianDOFBasis{D,T},
  prebasis::Basis{D,T}, b::AbstractMatrix{Float64}) where {D,T}
  vals = evaluate(prebasis,this.nodes)
  l = length(prebasis); lt = _length(T)
  # E = eltype(T)
  # b = Array{E,2}(undef,l, l)
  nnd = length(this.nodes)
  @assert nnd*_length(T) == length(prebasis)
  function computeb!(a,b,lt,nnd)
    for k in 1:lt
      off = nnd*(k-1)
      for j in 1:size(a,2)
        for i in 1:size(a,1)
          b[i,j+off] = a[i,j][k]
        end
      end
    end
  end
  computeb!(vals,b,lt,nnd)
  return b
end

"""
Evaluate the Lagrangian DOFs basis (i.e., nodal values) for a given field in
the reference space
"""
function numlocaldofs(this::LagrangianDOFBasis{D,T}) where {D,T}
  lt = _length(T)
  # E = eltype(T)
  nnd = length(this.nodes)
  # return b = Vector{E}(undef,lt*nnd)
  return lt*nnd
end

# @santiagobadia : Be careful, a physical field must be composed with geomap
# before being used here. Is this what we want?
function evaluate!(this::LagrangianDOFBasis{D,T},
  field::Field{D,T}, b::AbstractVector{Float64}) where {D,T}
  vals = Maps.evaluate(field,this.nodes)
  # I would like to use evaluate everywhere, putting evaluate in Gridap and
  # importing it in all submodules
  # This way we could use the same evaluate for bases and fields...
  # @santiagobadia : TO BE DONE
  lt = _length(T)
  # E = eltype(T)
  nnd = length(this.nodes)
  # b = Vector{E}(undef,lt*nnd)
  function computeb!(a,b,lt,nnd)
    for k in 1:lt
      off = nnd*(k-1)
      for j in 1:nnd
        b[j+off] = a[j][k]
      end
    end
  end
  computeb!(vals,b,lt,nnd)
  return b
end

_length(::Type{<:Real}) = 1

_length(::Type{T}) where T = length(T)

end # module
