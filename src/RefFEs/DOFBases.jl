module DOFBases

using Test
using Gridap
using Gridap.Helpers

export DOFBasis
export test_dof_basis
export numlocaldofs
export dof_type
export LagrangianDOFBasis

import Gridap: evaluate
import Gridap: evaluate!
import Base: length

# Interfaces

"""
Abstract DOF basis
"""
abstract type DOFBasis{D,T} end

"""
Returns the dimension of the DOF basis.
"""
function length(::DOFBasis)::Int
  @abstractmethod
end

"""
Compute the DOF values of the given field. The results in written in place
in the last argument.
E is a subtype of Real (typically Float64). E can be computed as
E = dof_type(T)
"""
function evaluate!(
  ::DOFBasis{D,T},::Field{D,T},::AbstractVector{E}) where {D,T,E}
  @abstractmethod
end

"""
Compute the DOF values of the given basis. The results in written in place
in the last argument.
E is a subtype of Real (typically Float64). E can be computed as
E = dof_type(T). The result is a matrix. The rows correspond to
the Basis and the cols to the DOFBasis.
"""
function evaluate!(
  ::DOFBasis{D,T},::Basis{D,T},::AbstractMatrix{E}) where {D,T,E}
  @abstractmethod
end

"""
Returns the type representing a DOF value. Typically Float64,
but can be other types depending on the input T.
"""
function dof_type(::Type{T}) where T
  E = eltype(T)
  @notimplementedif !(E<:Real)
  E
end

"""
Returns the dimension of the DOF basis.
"""
numlocaldofs(b::DOFBasis) = length(b)

function evaluate(b::DOFBasis{D,T},f::Field{D,T}) where {D,T}
  E = dof_type(T)
  n = length(b)
  v = zeros(E,n)
  evaluate!(b,f,v)
  v
end

function evaluate(b::DOFBasis{D,T},f::Basis{D,T}) where {D,T}
  E = dof_type(T)
  n = length(b)
  m = length(f)
  @assert n == m
  v = zeros(E,n,m)
  evaluate!(b,f,v)
  v
end

# Testers

function test_dof_basis(
  dofbasis::DOFBasis{D,T},
  field::Field{D,T},
  basis::Basis{D,T},
  field_dofs::AbstractVector{E},
  basis_dofs::AbstractMatrix{E}) where {D,T,E}

  @test dof_type(T) == E
  @test length(dofbasis) == length(basis)
  @test numlocaldofs(dofbasis) == length(basis)

  a = evaluate(dofbasis,field)
  @test a ≈ field_dofs

  b = evaluate(dofbasis,basis)
  @test b ≈ basis_dofs

end

# Concrete implementations

"""
Lagrangian DOF basis, which consists on evaluating the polynomial basis
(prebasis) in a set of points (nodes)
"""
struct LagrangianDOFBasis{D,T} <: DOFBasis{D,T}
  nodes::Vector{Point{D,Float64}}
  dof_to_node::Vector{Int}
  dof_to_comp::Vector{Int}
  node_and_comp_to_dof::Matrix{Int}
  _cache_field::Vector{T}
  _cache_basis::Matrix{T}
end

function LagrangianDOFBasis(::Type{T}, nodes::Vector{Point{D,Float64}}) where {D,T}

  nnodes = length(nodes)
  ndofs = _num_dofs(T,nodes)
  dof_to_node, dof_to_comp, node_and_comp_to_dof =
    _generate_dof_layout(T,nnodes)
  cache_field = zeros(T,nnodes)
  cache_basis = zeros(T,ndofs,nnodes)

  LagrangianDOFBasis{D,T}(
    nodes,
    dof_to_node,
    dof_to_comp,
    node_and_comp_to_dof,
    cache_field,
    cache_basis)
end

function length(b::LagrangianDOFBasis{D,T}) where {D,T}
  _num_dofs(T,b.nodes)
end

function evaluate!(
  b::LagrangianDOFBasis{D,T},f::Field{D,T},dofs::AbstractVector{E}) where {D,T,E}
  cache = b._cache_field
  evaluate!(f,b.nodes,cache)
  for (node,v) in enumerate(cache)
    for (comp,vi) in enumerate(v)
      dof = b.node_and_comp_to_dof[node,comp]
      dofs[dof] = vi
    end
  end
end

function evaluate!(
  b::LagrangianDOFBasis{D,T},f::Basis{D,T},dofs::AbstractMatrix{E}) where {D,T,E}
  cache = b._cache_basis
  evaluate!(f,b.nodes,cache)
  for node in 1:size(cache,2)
    for i in 1:size(cache,1)
      v = cache[i,node]
      for (comp,vk) in enumerate(v)
        dof = b.node_and_comp_to_dof[node,comp]
        dofs[i,dof] = vk
      end
    end
  end
end

# Helpers

# Node major implementation
function _generate_dof_layout(::Type{T},nnodes::Integer) where T
  ncomps = _length(T)
  ndofs = ncomps*nnodes
  dof_to_comp = zeros(Int,ndofs)
  dof_to_node = zeros(Int,ndofs)
  node_and_comp_to_dof = zeros(Int,nnodes,ncomps)
  v = zero(T)
  for node in 1:nnodes
    for comp in 1:ncomps
      o = nnodes*(comp-1)
      dof = node+o
      dof_to_comp[dof] = comp
      dof_to_node[dof] = node
      node_and_comp_to_dof[node,comp] = dof
    end
  end
  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

function _num_dofs(::Type{T}, nodes::Vector) where T
  ncomps = _length(T)
  nnodes = length(nodes)
  ncomps*nnodes
end

_length(::Type{<:Real}) = 1

_length(::Type{T}) where T = length(T)

end # module
