abstract type AbstractSymMultiValue{S,T,N,L} <: MultiValue{S,T,N,L} end

struct SymHighOrderTensorValue{S,T,N,L} <: AbstractSymMultiValue{S,T,N,L}
  data::NTuple{L,T}
end


