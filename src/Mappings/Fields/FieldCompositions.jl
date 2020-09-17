# Composition

Base.:∘(f::Field,g::Field) = GenericField(Composition(f,g))

struct Composition{F,G} <: Mapping
  f::F
  g::F
end

function return_cache(a::Composition,x)
  cg = return_cache(a.g,x)
  gx = testitem(a.g,x)
  cf = return_cache(a.f,gx)
  (cf, cg)
end

function evaluate!(cache,a::Composition,x)
  cf, cg = cache
  gx = evaluate!(cg,f.g,x)
  evaluate!(cf,f.f,gx)
end

function gradient(a::Composition)
  ∇f = gradient(a.f)
  ∇g = gradient(a.g)
  evaluate(FieldOperation(⋅),∇g,∇f)
end

# Composition of FieldArray and Field

Base.:∘(f::FieldArray,g::Field) = GenericFieldArray(Composition(f,g))
Base.size(a::Composition) = size(a.f)
Base.getindex(a::Composition,i::Integer...) = GenericField(Composition(a.f[i...],a.g))
Base.ndims(a::Composition{F}) where F = ndims(F)
Base.eltype(a::Composition{F,G}) where F = GenericField{Composition{eltype(F),G}}
Base.IndexStyle(::Type{<:Composition{F}}) where F = IndexStyle(F)
